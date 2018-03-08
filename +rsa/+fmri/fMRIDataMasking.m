function [varargout] = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence, userOptions);

% fMRIDataMasking is a function which takes a full set of fMRI brain scans and
% masks them according the the masks specified by the user (if any).  If no
% masks are specified then no masking is done.
%
% [responsePatterns =] fMRIDataMasking(fullBrainVols, binaryMasks_nS, userOptions)
%
%       fullBrainVols --- The unmasked beta (or t) images.
%               fullBrainVols.(subject) is a [nVoxel nCondition nSession]-sized
%               matrix. The order of the voxels is that given by reshape or (:).
%
%        binaryMasks_nS --- The native- (subject-) space masks.
%               binaryMasks_nS.(subject).(mask) is a [x y z]-sized binary matrix
%               (the same size as the native-space 3D beta images.
%
%        betaCorrespondence --- The array of beta filenames.
%               betas(condition, session).identifier is a string which referrs
%               to the filename (not including path) of the SPM beta image. (Or,
%               if not using SPM, just something, as it's used to determine the
%               number of conditions and sessions.)
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
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Defaults to the fieldnames of the first subject
%                        of binaryMasks_nS.
%
% The following files are saved by this function:
%        userOptions.rootPath/ImageData/
%                userOptions.analysisName_responsePatterns.mat
%                        Contains the masked brains in a struct such that
%                        responsePatterns.(mask).(subject) is a [nMaskedVoxels
%                        nConditions nSessions] matrix.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_fMRIDataMasking_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
%  Cai Wingfield 12-2009, 6-2010

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('fMRIDataMasking:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('fMRIDataMasking:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(fullBrainVols));
userOptions = setIfUnset(userOptions, 'maskNames', fieldnames(binaryMasks_nS.(userOptions.subjectNames{1})));

% The analysisName will be used to label the files which are eventually saved.
MaskedBrainsFilename = [userOptions.analysisName, '_responsePatterns.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRIDataMasking_Details.mat'];

promptOptions.functionCaller = 'fMRIDataMasking';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag

	% Data
	nMasks = numel(userOptions.maskNames);

	%% Get Data

	if ischar(betaCorrespondence) && strcmpi(betaCorrespondence, 'SPM')
        betas = getDataFromSPM(userOptions);
    else
        betas = betaCorrespondence;
    end%if:SPM

	nSubjects = numel(userOptions.subjectNames);
	nSessions = size(betas, 1);
	nConditions = size(betas, 2);

	for mask = 1:nMasks % For each mask...

		% Figure out which mask we're looking at

		thisMask = userOptions.maskNames{mask};

		nanFlagOut = zeros(nSubjects); % This is a matrix which will store any time when a subject's masks peek outside the brain, and if they're inconsistent between scans (in which case results are called into some question).

		for subject = 1:nSubjects % and for each subject...

			% Figure out which subject this is
			thisSubject = userOptions.subjectNames{subject};

			% Get the brain scan vol
			thisScan = fullBrainVols.(thisSubject);

			% Load the mask data into a vector
			maskMatrix = binaryMasks_nS.(thisSubject).(thisMask);

			maskVector = reshape(maskMatrix, 1, []);

			% Check to see if the mask is empty
			if numel(find(maskVector)) == 0
				error(['The mask "' thisMask '" for subject "' thisSubject '" is empty!  Check your mask definitions.']);
			end%if

			nanFlagIn = zeros(nSessions, nConditions); % This is a matrix which, for each subject will store any time when NaNs are removed from subject's RoI between sessions or conditions.

			for session = 1:nSessions % and each session...

				% Apply the mask
				whereToLook = find(maskVector); % Inside RoI
				for voxel = 1:max(size(whereToLook))
					insideMaskVoxelsWithNans(voxel, :, session) = thisScan(whereToLook(voxel), :, session);
				end%for

				%% Now avoid NaNs introduced by the reverse-normalisation of anatomically-defined RoI masks (or full brains!?)
				for condition = 1:nConditions

					% nanFlagIn provides information about whether NaNs have been removed, and if it's more than once.
					if numel(find(isnan(insideMaskVoxelsWithNans(:, condition, session)))) > 0
						nanFlagIn(session, condition) = nanFlagIn(session, condition) + 1;
					end%if
					% By the end of this, if NaNs have been removed once, nanFlagIn should be 1, and if they've been removed more than once, nanFlagIn should be 2 or more.

					whereToLook = find(~isnan(insideMaskVoxelsWithNans(:, condition, session)));

					for voxel = 1:max(size(whereToLook))
						insideMaskVoxels(voxel, condition, session) = insideMaskVoxelsWithNans(whereToLook(voxel), condition, session);
					end%for
				end%for
			end%for

			% Put the masked data into a struct (to be saved)
			responsePatterns.(thisMask).(thisSubject) = insideMaskVoxels;

			clear insideMaskVoxels insideMaskVoxelsWithNans thisScan whereToLook;

			%% Warn user about inconsistently masked data if necessary

			nanFlagOut(subject) = sum(sum(nanFlagIn)); % This should now store the number of times when NaNs had to be removed

			if nanFlagOut(subject) > 0
				fprintf([ ...
					'+\n' ...
					'+\tThe mask "' thisMask '" caught some NaNs in subject\n' ...
					'+\t"' thisSubject '", probably looking outside the brain.\n' ...
					'+\tResults are therefore not based on all of\n' ...
					'+\t"' thisMask '".\n']);
			end%if

			if nanFlagOut(subject) > 1
				fprintf([ ...
					'+\n' ...
					'+\tNote: This is not the first time it has happened with\n' ...
					'+\t      this mask [x' num2str(nanFlagOut(subject)) '].  Brains or masks may differ\n' ...
					'+\t      between runs. Maybe check your data?  Interpret\n' ...
					'+\t      results for "' thisSubject '" and "' thisMask '"\n' ...
					'+\t      at your own risk!\n']);
			end%if

		end%for
	end%for

	%% Save relevant info

	timeStamp = datestr(now);

	fprintf(['Saving masked data to ' fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'ImageData');
	save(MaskedBrainsFilename, 'responsePatterns','-v7.3');

	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');
	
else
	fprintf(['Loading previously saved RoIs from ' fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'ImageData', MaskedBrainsFilename));
end%if

if nargout == 1
	varargout{1} = responsePatterns;
elseif nargout > 0
	error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started

end%function
