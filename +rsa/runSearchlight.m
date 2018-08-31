function runSearchlight(Searchlight,inFiles,outFiles,rsaFunction,varargin)
% function runSearchlight(Searchlight,infiles,outfiles,rsaFunction,varargin)
% Main analysis engine for searchlight-based analysis
% The function deals with data loading from the infiles, concatinates the
% data to a nObservation x nVoxel matrix, and then hands these on to the
% rsaFunction 
% 
% INPUT:
%   Seachlight:  Search light definition. This is a structure that contains 
%               LI:         nVox x 1 cell array with linear voxel indices
%               voxmin:     nVox x 3 Minimal voxelcoordinate in x,y,z direction
%               voxmax:     nVox x 3 maximal voxelcoordinate in x,y,z direction
%               voxel:      nVox x 3 Matrix of I,J,K voxel coordinates or 
%                           nVox x 1 vector of linear voxel indices for the centers of search-lights 
%   infiles:     Input images file names or memory-mapped volumes (SPM) 
%   outfiles:    Output file name (in cell array)
%   rsaFunction: The multivariate analysis function, it should take X (data, N x P)
%                as the first input argument, and then receives any other
%                optionalParams. The first output of the function should be a Hx1 or 1xH
%                vector, which will be stored in the output images 
% USEROPTIONS / VARARGIN:
%   optionalParams: extra parameters passed to the mva-function
%   subset:         So the analysis only on subset of voxels 
%   transposeInput: allows the rsaFunction to take a nVox x nObs matrix 
%   idealBlock:     ideal Block size (default 7e5) 
% (C) Joern Diedrichsen 2015 

import rsa.surf.*

% User optional parameters 
Opt.optionalParams = []; 
Opt.subset         = []; 
Opt.transposeInput = 0; 
Opt.idealBlock     = 7e5;    % Ideal block size, increase this if you have more free memory

Opt = rsa.getUserOptions(varargin,Opt); 

% Check input files 
if (iscell(inFiles))
    inFiles=char(inFiles);
end;
if (isstruct(inFiles))
    VolIn=inFiles; 
else 
    VolIn=spm_vol(inFiles);
end; 

% Check the searchlight input
if (~isstruct(Searchlight) || ~isfield(Searchlight,'LI') || ~isfield(Searchlight,'voxel')...
                           || ~isfield(Searchlight,'voxmin') || ~isfield(Searchlight,'voxmax'));
    error('searchlight structure needs to be a Structure with fields .LI, .vox, .voxmin, .voxmax'); 
end;
if (~iscell(Searchlight.LI) || size(Searchlight.LI,2)~=1)
    error('Searchlight.LI must be a nVox x 1 cell array');
end;

% determine length of output voxels 
numOutfiles = length(outFiles);
numVoxel    = size(Searchlight.LI,1);
numInfiles  = length(inFiles);
% Check output files 
if (~iscell(outFiles))
    error('outfiles must be given in a cell array'); 
end; 

% Check if rsaFunction gives the correct number of outputs
X = rand(numInfiles,100); 
if (Opt.transposeInput) 
    Result = feval(rsaFunction,X',Opt.optionalParams{:});
else 
    Result = feval(rsaFunction,X,Opt.optionalParams{:});
end; 
if (length(Result))~=numOutfiles
    error('1st return argument of the rsaFunction must be a vector with the same length as the number of output files'); 
end; 

% Check subset option 
if (isempty(Opt.subset))
    Opt.subset = true(numVoxel,1); 
end; 

% Check 
if (size(Searchlight.voxel,2)==3)       % I,J,K coordinates 
    Searchlight.voxel = rsa.surf.surfing_subs2inds(VolIn(1).dim,Searchlight.voxel(:,1:3));
elseif (size(Searchlight.voxel,2)==1)   % Linear coordinates: do nothing 
else 
    error('Searchlight.voxel needs to be nVoxel x 3 or nVoxel x 1 ');
end; 

% Now split the search volume in multiple Blocks and save resulting
% structures as temp files to preserve memory
isDone      =  false;
isCalc      =  isnan(Searchlight.voxmin(:,1)) | ~Opt.subset;      % Set all NaN to calculated
blockSize   =  max(Searchlight.voxmax-Searchlight.voxmin);        % Minimal block size 
ratio       =  Opt.idealBlock./(prod(blockSize)*length(inFiles)); % See if the idealBlock size is bigger 
if (ratio>1)
   blockSize = blockSize*ratio^(1/3);                             % Increase blocksize in all directions 
end; 
b=1;

while (~isDone)
    % Find the candidate search lights
    % This starts in the x-direction and finds the slice of search lights
    % that are within that slice, and similar for y and z.
    % then it moves the block tightly to this location
    isCandidate = ~isCalc;
    for i=1:3
        IJKc1(i)    = min(Searchlight.voxmin(isCandidate,i));
        IJKc2(i)    = IJKc1(i)+blockSize(i);     
        isCandidate = isCandidate & Searchlight.voxmin(:,i)>=IJKc1(i) & ...
                                    Searchlight.voxmax(:,i)<=IJKc2(i);
    end;
    IJKc1  = min(Searchlight.voxmin(isCandidate,:),[],1);
    IJKc2  = IJKc1+blockSize;
    isCandidate=~isCalc & Searchlight.voxmin(:,1)>=IJKc1(1) & Searchlight.voxmin(:,2)>=IJKc1(2) &  ...
                          Searchlight.voxmin(:,3)>=IJKc1(3) & Searchlight.voxmax(:,1)<=IJKc2(1) &  ...
                          Searchlight.voxmax(:,2)<=IJKc2(2) & Searchlight.voxmax(:,3)<=IJKc2(3);
    j=find(isCandidate);
    if (isempty(j));
        break;
    end;
    
    % Save the substructure as a tempory file
    T.LI={Searchlight.LI{j}}';
    T.j=j;
    T.voxel=Searchlight.voxel(j,:);
    fprintf('block %d Corner: %d %d %d length:%d  \n',b,IJKc1,length(j));
    save(sprintf('temp_%2.2d.mat',b),'T');
    isCalc(j) = true;
    isDone    = all(isCalc);
    b=b+1;
end
numBlocks=b-1;

% Free Some memory
clear S T;

% Initialize output arrays
Result=zeros(numOutfiles,numVoxel)*NaN;

% Initialize progress bar and start timer 
k=1;
spm_progress_bar('Init',numBlocks,'MVA','number of Blocks');
tic;

% Now iterate over the number of blocks 
for b=1:numBlocks;
    load(sprintf('temp_%2.2d.mat',b));
    
    % Find the linear voxel indices 
    linVox  = unique(cat(2,T.LI{:})');
    [I,J,K] = ind2sub(VolIn(1).dim,linVox);
    
    % Get the data 
    X = sparse(double(max(linVox)),length(VolIn));
    N = sparse(double(max(linVox)),1);
    for i=1:length(VolIn)
        X(linVox,i)=spm_sample_vol(VolIn(i),double(I),double(J),double(K),0);
    end;
    clear I J K; % Keep memory small 
    
    % Now call the rsaFunction with the input 
    for i = 1:size(T.LI,1);
        if (Opt.transposeInput) 
            Result(:,T.j(i)) = feval(rsaFunction,full(X(T.LI{i},:)),Opt.optionalParams{:});
        else 
            Result(:,T.j(i)) = feval(rsaFunction,full(X(T.LI{i},:)'),Opt.optionalParams{:});
        end; 
        if (mod(k,50)==0)
            fprintf('%d done: %f\n',k,toc);
            tic;
        end;
        k=k+1;
    end;
    spm_progress_bar('Set',b);
    b=b+1;
end;
fprintf('%d total: %f\n',k,toc);

% Clean up the search lights
delete('temp*.mat');


% Ensure that the descriptor is a string 
if (~ischar(rsaFunction))
    rsaFunction = func2str(rsaFunction); 
end; 

% Now write output as a volume-based image 
for i=1:size(Result,1)
    Z=NaN(VolIn(1).dim);
        [d f t m]=spm_fileparts(outFiles{i});
    
    Vo      = struct(...
            'fname',    fullfile(d,[f t]),...
            'dim',      VolIn(1).dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      VolIn(1).mat,...
            'n',        [str2num(m) 1],...
            'descrip',  rsaFunction);
    Z(Searchlight.voxel)=Result(i,:);            % Project back into voxel space 
    spm_write_vol(Vo,Z);
end;
spm_progress_bar('Clear');

