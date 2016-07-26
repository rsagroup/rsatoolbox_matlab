function runSearchlightLDC(searchLight,varargin)
% runSearchlightLDC(searchLight,varargin)
% Wrapper for the main searchlight function rsa_runSearchlight 
% This version calculates the LDC from an SPM first-level analysis 
% New Version takes into account the first-level design matrix
% 
% INPUT: 
%     searchLight:  a) Name of the precalculated searchlight file 
%                   b) Searchlight structure itself (see rsa_defineSearchlight.m) 
% 
% VARARGIN: either list of 'optionname',value,'optioname2',value2.. 
%           or structure with Opt.optioname = value 
%   rootPath:       Path to be pre-pended to all path 
%   spmDir:         Path (relative to root) where the SPM.mat resides and
%                   were results are written 
%   conditionVec:   Vector of condition labels for the task-related
%                   regressors in the SPM matrix, if empty, the routine
%                   assumes that all betas are condition1.. conditionK
%   conditionLabels:Cell array of condition labels. If given, the routine
%                   attempts to extract the conditions and partitions from
%                   the SPM structure 
%   imageDataDir:   Alternative directory to find the raw imaging data, if
%                   it different from the one that is specified in the SPM
%                   structure. This is especially useful if the SPM was
%                   estimated in one directory and then later moved. 
% (C) Joern Diedrichsen 2015 

import rsa.spm.*
import rsa.util.*

% User optional parameters 
Opt.rootPath        = []; 
Opt.spmDir          = [];     %  Directory that contains the first-level analysis 
Opt.analysisName    = [];     %  Name of the analysis 
Opt.conditionVec    = [];     %  Vector, indicating which of the betas belong to which condition 
                              %  Use 0 for betas that will not be included 
Opt.conditionLabels = {};     %  Condition labels - used if conditionVec is not given
Opt.partition       = [];     %  Indicator of partition 
Opt.imageDataDir      = [];   %  Directory where the pre-processed raw data resides (if different from what is specified in the SPM-structure) 
Opt.saveSigma       = 'none'; %  Determine whether to save the variance / covariance of the beta estimates as well as the distances 

Opt = rsa.getUserOptions(varargin,Opt); 

workingDir = fullfile(Opt.rootPath, Opt.spmDir); 

% get the searchlight 
if (isstruct(searchLight)) 
elseif (ischar(searchLight)) 
    [dir,name,ext] = fileparts(searchLight);
    if (isempty(dir))
        searchLight = fullfile(workingDir,searchLight); 
    end; 
    searchLight = load(searchLight); 
end; 

% Load the SPM structure 
load(fullfile(workingDir,'SPM.mat')); 

% Determine the condition and the partition labels 
nBetas = length(SPM.xX.iC);                         % Use only task related regressors (no intercepts) 
if (isempty(Opt.conditionVec))
    if (isempty(Opt.conditionLabels))
        error('Either Opt.conditionVec or Opt.conditionLabels needs to be set'); 
    end;
    [Opt.conditionVec,Opt.partition]=rsa.getSPMconditionVec(SPM,Opt.conditionLabels); 
else 
    if (isempty(Opt.partition))
        Opt.partition=zeros(nBetas,1); 
        for i=1:length(SPM.Sess)
            Opt.partition(SPM.Sess(i).col)=i; 
        end; 
    end; 
end; 

nConditions = max(Opt.conditionVec); 

% Now define the input files 
if (isempty(Opt.imageDataDir)) 
    inFiles = SPM.xY.VY; 
else 
    numFiles = length(SPM.xY.VY);
    for i=1:numFiles
        [dir,filename,extension,number]=spm_fileparts(SPM.xY.P(i,:));
        inFileNames{i} = fullfile(Opt.rootPath,Opt.imageDataDir,sprintf('%s%s%s',filename,extension,number));
    end; 
    inFiles = spm_vol(char(inFileNames)); 
end; 

% Now define the outfiles 
nDistances  =  nConditions * (nConditions-1)/2;
for i=1:nDistances 
    outFiles{i} = fullfile(workingDir,sprintf('%s_LDC.nii,%d',Opt.analysisName,i)); 
end; 

% Now generate the design condition matrix for the betas 
C   = rsa.util.indicatorMatrix('allpairs',[1:nConditions]);
Z   = rsa.util.indicatorMatrix('identity_p',Opt.conditionVec);

% Finally call the main searchlight engine 
rsa.runSearchlight(searchLight,inFiles,outFiles,@calculateLDC,'optionalParams',{SPM,Opt.partition,Opt.conditionVec});

% This is the plug-in function that is being called from rsa_runSearchlight 
function output = calculateLDC(Y,SPM,partition,conditionVec); 
U_hat   = rsa.spm.noiseNormalizeBeta(Y,SPM);          % Get noise normalised betas 
output  = rsa.distanceLDC(U_hat,partition,conditionVec,SPM.xX.xKXs.X);              % record distances as output 
