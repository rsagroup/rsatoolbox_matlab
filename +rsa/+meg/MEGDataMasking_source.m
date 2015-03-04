% MEGDataMasking_source is a module which masks source-localised EMEG
% time-course data to RoIs specified by lists of time- and vertex-indices.
%
% [maskedMeshes =] MEGDataMasking_source(
%                                        sourceMeshes,
%                                        indexMasks,
%                                        maskSpec, % updated to searchlightPatterns IZ
%                                        betaCorrespondence,
%                                        userOptions
%                                       )
%
%        sourceMeshes --- The unmasked sourceMeshes topographies.
%               Contains the subject source-reconstructed mesh data in a
%               struct such that sourceMeshes.(subjectName).(L|R) is a
%               [nVertices nTimepoints nConditions nSessions]-sized matrix.
%
%        indexMasks --- The specification of the masks which will be applied
%                       to the data.
%               indexMasks.(maskName).maskIndices
%                       A vector of indices inside the mask (indices above
%                       10242 are ignored).
%               indexMasks.(maskName).timeIndices
%                       A vector of indices for timepoints to be included
%                       within the mask.
%               indexMasks.(maskName).chirality
%                       Either "L" or "R", depending on which hemisphere the
%                       mask indices refer to. Cross-hemisphere RoIs aren't
%                       currently supported...
%   
%        maskSpec --- The specifications of the masking. % update IZ: see
%                       userOptions.searchlightPatterns for details
%               maskSpec.patternType
%                       What kind of patterns should be kept? Can take the
%                       following values:
%                               'Spatial'
%                                       A median is taken over the time window
%                                       of interest so RDMs will be calculated
%                                       based on sensor patterns over space.
%                               'Temporal'
%                                       All sensors are (mean) averaged over so
%                                       RDMs will be calculated based on
%                                       patterns over time.
%                               'Spatiotemporal'
%                                       No averaging is done, the time-courses
%                                       for each sensor with in the time window
%                                       of interest are concatenated. RDMs will
%                                       be calculated based on patterns over
%                                       both time and space.
%                       Defaults to 'Spatial'.
%
%        betaCorrespondence --- The array of beta filenames.
%                betas(condition, session).identifier is a string which referrs
%                to the filename (not including path) of the SPM beta image.
%                (Or, if not using SPM, just something, as it's used to
%                determine the number of conditions and sessions.)
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in sensorImages.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Defaults to the fieldnames of indexMasks.
%
% The following files are saved by this function:
%        userOptions.rootPath/ImageData/
%                userOptions.analysisName_MaskedSensors.mat
%                        Contains the raw timecourses in a structure such that
%                        maskedSensors.masked.(subjectName) is a [nChannels
%                        nConditions nSessions]-sized matrix (median across time
%                        window).
%        userOptions.rootPath/Details/
%                userOptions.analysisName_MEGDataMasking_sensor_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%  
% Cai Wingfield 9-2010
% Updated Isma Zulfiqar 11-2012
% Updated Fawad 26/03/2014

function [varargout] = MEGDataMasking_source(sourceMeshes, indexMasks, betaCorrespondence, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

maximumVertexIndex = userOptions.nVertices;

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('MEGDataMasking_source:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('MEGDataMasking_source:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(sourceMeshes));
userOptions = setIfUnset(userOptions, 'maskNames', fieldnames(indexMasks));
% maskSpec = setIfUnset(maskSpec, 'patternType', 'Spatial'); %removed due to new options

% The analysisName will be used to label the files which are eventually saved.
MaskedBrainsFilename = [userOptions.analysisName, '_MaskedMeshes.mat'];
DetailsFilename = [userOptions.analysisName, '_MEGDataMasking_source_Details.mat'];

promptOptions.functionCaller = 'MEGDataMasking_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

if ~userOptions.slidingTimeWindow, overwriteFlag = overwritePrompt(userOptions, promptOptions); 
else overwriteFlag = 1; end

if overwriteFlag

	% Data
	nMasks = numel(userOptions.maskNames);

	%% Get Data

	betas = betaCorrespondence;

	nSubjects = numel(userOptions.subjectNames);
	nSessions = size(betas, 1);
	nConditions = size(betas, 2);

	for mask = 1:nMasks
	
		% Which mask is this?
		thisMask = dashToUnderscores(userOptions.maskNames{mask}); % updated by Li Su 2-2012
        
          
            
            
            % Load the mask data into a vector
            maskIndices = indexMasks.(thisMask).maskIndices;
            maskIndices = maskIndices(maskIndices <= maximumVertexIndex);
            % timeIndices = indexMasks.(thisMask).(thisTimeWindow).timeIndices;
            timeIndices = userOptions.maskTimetoDataPoints.(thisMask); % updated IZ 11-12
            chi = indexMasks.(thisMask).chirality;

          for subject = 1:nSubjects % and for each subject...

                % Figure out which subject this is
                thisSubject = userOptions.subjectNames{subject};

                % Get the mesh for this subject
                thisMesh = sourceMeshes(subject).(chi);

                % Mask the data
                if userOptions.nSessions==1
                    maskedMesh = thisMesh(maskIndices, timeIndices(1):timeIndices(2), :); % (vertices, timePointes, conditions) % updated IZ 11-12
                                       
                else
                    maskedMesh = thisMesh(maskIndices, timeIndices(1):timeIndices(2), :, :); % (vertices, timePointes, conditions, sessions) % updated IZ 11-12
                    
                end
             if ~userOptions.regularized
                % Reduce to correct data type
                switch lower(userOptions.searchlightPatterns)
                    % For spatial patterns, median across the time window
                    case 'spatial'
                        reducedMaskedMesh = zeros(size(maskedMesh, 1), size(maskedMesh, 3), size(maskedMesh, 4)); % (data, conditions, sessions)
                        reducedMaskedMesh(:,:,:) = median(maskedMesh, 2);
                    % For temporal patterns, average across the vertices.
                    case 'temporal'
                        reducedMaskedMesh = zeros(size(maskedMesh, 2), size(maskedMesh, 3), size(maskedMesh, 4)); % (data, conditions, sessions)
                        reducedMaskedMesh(:,:,:) = mean(maskedMesh, 1);
                    % For spatiotemporal patterns, concatenate the vertex timecourses
                    case 'spatiotemporal'
%                         reducedMaskedMesh = zeros(1, size(maskedMesh, 3), size(maskedMesh, 4)); % (data, conditions, sessions)
%                         for vertex = 1:size(maskedMesh, 1)
%                             if vertex == 1
%                                 reducedMaskedMesh(:,:,:) = maskedMesh(vertex, :, :, :);
%                             else
%                                 reducedMaskedMesh = [reducedMaskedMesh; maskedMesh(vertex, :, :, :)];
%                             end%if
%                         end%for:vertex
                        reducedMaskedMesh = reshape(maskedMesh, [], size(maskedMesh, 3), size(maskedMesh, 4)); % update IZ 11-12
                end % switch
                
              else% update IZ 11-12 
                 % case 'regularized' 
                 tempMesh = reshape(maskedMesh, [], size(maskedMesh, 3), size(maskedMesh, 4));   
                 reducedMaskedMesh = zeros(size(maskedMesh, 1)*size(maskedMesh,2), size(maskedMesh, 3)* size(maskedMesh, 4)); % (data, conditions, sessions)

                  % combining session-wise trials
                 k=1;
                 for j=1:size(tempMesh,2)
                    for i=1:nSessions
                        reducedMaskedMesh(:,k) = (tempMesh(:,j,i)); 
                        k=k+1;
                    end
                 end
              end % if regularized                             

                % Store in struct
                maskedMeshes.(thisMask).(thisSubject) = reducedMaskedMesh;
                % maskedMeshes.(thisMask) = reducedMaskedMesh; % IZ 11/12

            end%for:subject
	end%for:mask

	%% Save relevant info

	timeStamp = datestr(now);
    
	fprintf(['Saving masked data to ' fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'ImageData');
	save(MaskedBrainsFilename, 'maskedMeshes');

	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');%, 'searchlightPatterns');
	
else
	fprintf(['Loading previously saved RoIs from ' fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename));
end%if

if nargout == 1
	varargout{1} = maskedMeshes;
elseif nargout > 0
	error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started
