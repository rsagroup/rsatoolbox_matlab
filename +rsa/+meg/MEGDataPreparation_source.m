% MEGDataPreparation_source is a function designed to take MEG data and load it
% into the correct format for the rest of the toolbox to use.
%
% It is based on Su Li's code
%
% [sourceMeshes[, baselineLimit] =] MEGDataPreparation_source(
%                                                            subject,
%                                                            betaCorrespondence,
%                                                            userOptions)
%
%       betaCorrespondence --- The array of beta filenames.
%               betas(condition, session).identifier is a string which referrs
%               to the filename (not including path) of the SPM beta image.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in fullBrainVols.
%                userOptions.betaPath
%                        A string which contains the absolute path to the
%                        location of the beta images. It can contain the
%                        following wildcards which would be replaced as
%                        indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%                                [[betaIdentifier]]
%                                        To be replaced by filenames as provided
%                                        by betaCorrespondence.
%
% The following files are saved by this function:
%        userOptions.rootPath/ImageData/
%                userOptions.analysisName_CorticalMeshes.mat
%                        Contains the subject source-reconstructed mesh data in
%                        a struct such that sourceMeshes.(subjectName).(L|R) is
%                        a [nVertices nTimepoints nConditions nSessions]-sized
%                        matrix.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_MEGDataPreparation_source_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
% REQUIRES MATLAB VERSION 7.3 OR LATER!
%
% Li Su's Note on 2-2012: the previous version of this function is not time
% and space effective because when there are too many conditions, the
% sourceMashes is too big to fit in the memory, and it takes a lot of disk
% space. And the saved sourceMashes is slow when loading into the memory.
% So, the current version will constract sourceMashes for each subject on
% the fly and only store the r-map for each model fitting instead of saving
% sourceMashes on to the disk. When the condition number is small, it
% better to save the data RDMs to the disk so that we don't need to
% recompute it when comparing it to differnt model. But for most of the
% time, it is not possible because loading them from the disk is still
% slow. As the reasons above, I decided to change this function from a
% module to an engine function, and call it from MEGSearchlight_source.m 
% and make it subject specific too. 
%
%  Cai Wingfield 5-2010, 6-2010 updated by Li Su 2-2012 
%   updated Fawad 3-12014


function [varargout] = MEGDataPreparation_source(subject, betaCorrespondence, userOptions)

% returnHere = pwd; % We'll return to the pwd when the function has finished
% 
% %% Set defaults and check options struct
% 
% % The filenames contain the analysisName as specified in the user options file
% DetailsFilename = [userOptions.analysisName, '_MEGDataPreparation_source_Details.mat'];
% BaselineLimitFilename = [userOptions.analysisName, '_BaselineLimit.mat'];
% 
% 
% promptOptions.functionCaller = 'MEGDataPreparation_source';
% promptOptions.defaultResponse = 'S';
% promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);
% ImageDataFilename = [userOptions.analysisName, '_', userOptions.subjectNames{nSubjects}, '_CorticalMeshes.mat']; % only check the LAST subject
% promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename);
% 
% overwriteFlag = overwritePrompt(userOptions, promptOptions);
% 
% if overwriteFlag

	%% Get Data
	betas = userOptions.betaCorrespondence;
	[nSessions, nConditions] = size(betas);
    downsampleRate = userOptions.temporalDownsampleRate;
    
    % First get just the first condition of the first subject's left data,
    % to set the correct sizes, even if a future reading fails due to
    % artifacts
    readPathL = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', betas(1, 1).identifier, '[[subjectName]]', userOptions.subjectNames{1}, '[[LR]]', 'l');
    MEGDataStcL = mne_read_stc_file1(readPathL);
    MEGDataStcL.data = MEGDataStcL.data(:,1:downsampleRate:end);
    [nVertex_Raw, nTimepoint_Raw] = size(MEGDataStcL.data);
    
	%for subject = 1:nSubjects % For each subject

		% Figure out the subject's name
		thisSubject = userOptions.subjectNames{subject};
        missingFilesLog = fullfile(userOptions.rootPath, 'ImageData', 'missingFilesLog.txt');
        if ~exist(fullfile(userOptions.rootPath, 'ImageData'),'dir')
            mkdir(userOptions.rootPath,'ImageData')
        end
        %ImageDataFilename = [userOptions.analysisName, '_', thisSubject, '_CorticalMeshes.mat'];
		
		for session = 1:nSessions % For each session...  UNUSED?!?!
			for condition = 1:nConditions % and each condition...

				% Then read the brain data (for this session, condition)
				readPathL = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', betas(session, condition).identifier, '[[subjectName]]', thisSubject, '[[LR]]', 'l');
				readPathR = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', betas(session, condition).identifier, '[[subjectName]]', thisSubject, '[[LR]]', 'r');
				try
                    MEGDataStcL = mne_read_stc_file1(readPathL);
                    MEGDataStcR = mne_read_stc_file1(readPathR);
                    % data reordering by vertex list IZ 06/13
                    MEGDataStcL.data = orderDatabyVertices(MEGDataStcL.data, MEGDataStcL.vertices);
                    MEGDataStcR.data = orderDatabyVertices(MEGDataStcR.data, MEGDataStcR.vertices);
                    MEGDataVolL = single(MEGDataStcL.data);
                    MEGDataVolR = single(MEGDataStcR.data);
                    MEGDataVolL = MEGDataVolL(:,1:downsampleRate:end);
                    MEGDataVolR = MEGDataVolR(:,1:downsampleRate:end);
                catch ex
                    % when a trial is rejected due to artifact, this item
                    % is replaced by NaNs. Li Su 3-2012
                    disp(['Warning: Failed to read data for condition ' num2str(condition) '... Writing NaNs instead.']);
                    dlmwrite(missingFilesLog, str2mat(replaceWildcards( betas(session, condition).identifier, '[[subjectName]]', thisSubject)), 'delimiter', '', '-append');
                    MEGDataVolL = NaN(nVertex_Raw,nTimepoint_Raw);
                    MEGDataVolR = MEGDataVolL;
				end
				
				baselineLimit = double(-MEGDataStcL.tmin/MEGDataStcL.tstep); % I hope these are all the same!
                
                subjectSourceData.L(:, :, condition, session) = MEGDataVolL; % (vertices, time, condition, session)
				subjectSourceData.R(:, :, condition, session) = MEGDataVolR;
                
                if nConditions > 10
                    if mod(condition, floor(nConditions/20)) == 0, fprintf('\b.:'); end%if
                else
                    fprintf('\b.:');
                end
				
			end%for
            % fprintf('\b.:');
            
		end%for
        fprintf('\bData read successfully!\n');
        dlmwrite(missingFilesLog, '', '-append');

		% For each subject, record the vectorised brain scan in a subject-name-indexed structure
		
		%% MEMORY DEBUG %%
		sourceMeshes = subjectSourceData;
%  		cd(fullfile(userOptions.rootPath, 'ImageData'));
%         fprintf(['Saving image data to ' fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename) '\n']);
%  		save(ImageDataFilename, 'sourceMeshes'); 		
%		clear sourceMeshes subjectSourceData MEGDataStcL MEGDataStcR MEGDataVolL MEGDataVolR;

	%end%for:subjects

	%% Save relevant info

% 	timeStamp = datestr(now);
% 	%save(ImageDataFilename, 'sourceMeshes', '-v7.3'); % -v7.3 required as this is likely to be larger than 2GB
% 	%% END DEBUG %%
% 	save(BaselineLimitFilename, 'baselineLimit');
% 	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
% 	gotoDir(userOptions.rootPath, 'Details');
% 	save(DetailsFilename, 'timeStamp', 'userOptions');
	
% else
%     fprintf(['Loading previously saved meshes from ' fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename) '...\n']);
%     load(fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename)); 
% 	fprintf(['Loading previously saved baseline limit from ' fullfile(userOptions.rootPath, 'ImageData', BaselineLimitFilename) '...\n']);
% 	load(fullfile(userOptions.rootPath, 'ImageData', BaselineLimitFilename));
%end%if

if nargout == 1
	varargout{1} = sourceMeshes;
elseif nargout == 2
	varargout{1} = sourceMeshes;
	varargout{2} = baselineLimit;
elseif nargout > 0
	error('0, 1 or 2 arguments out, please.');
end%if:nargout

%cd(returnHere); % Go back
