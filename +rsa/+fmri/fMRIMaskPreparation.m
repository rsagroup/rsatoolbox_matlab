function [varargout] = fMRIMaskPreparation(userOptions)
%
% fMRIMaskPreparation will load SPM-defined masks from a directory specified in
% userOptions and will save a struct containing binary mask matrices.
%
% [binaryMasks_nS =] fMRIMaskPreparation(userOptions)
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.maskPath
%                        A string describing the path to the location of the
%                        files for the definition of RoI masks.
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names.
%
% The following files are saved by this function:
%        userOptions.rootPath/ImageData/
%                userOptions.analysisName_Masks.mat
%                        Contains the raw beta images in a struct such that
%                        binaryMasks_nS.(subject).(mask) is a [x y z]-sized
%                        binary matrix.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_fMRIMaskPreparation_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
% Cai Wingfield 2-2010, 6-2010

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
if ~isfield(userOptions, 'analysisName'), error('fMRIMaskPreparation:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('fMRIMaskPreparation:NoRootPath', 'rootPath must be set. See help'); end%if
if ~isfield(userOptions, 'subjectNames'), error('fMRIMaskPreparation:NoSubjectNames', 'betaPath must be set. See help'); end%if
if ~isfield(userOptions, 'maskNames'), error('fMRIMaskPreparation:NoMaskNames', 'maskNames must be set. See help'); end%if
if ~isfield(userOptions, 'maskPath'), error('fMRIMaskPreparation:NoMaskPath', 'maskPath must be set. See help'); end%if

% The analysisName will be used to label the files which are eventually saved.
MasksFilename = [userOptions.analysisName, '_Masks.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRIMaskPreparation_Details.mat'];

promptOptions.functionCaller = 'fMRIMaskPreparation';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'ImageData', MasksFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag

	% Data
	nSubjects = numel(userOptions.subjectNames);
	nMasks = numel(userOptions.maskNames);
	
	% Check to make sure it's not empty
	%if isempty(masks), error('(!) No masks were selected in userOptions.maskChoices\n'); end%if
	
	%% Get Data
	
	fprintf('Loading binary masks...\n');

	for subjectNumber = 1:nSubjects % and for each subject...

		% Figure out which subject this is
		subject = userOptions.subjectNames{subjectNumber};
		
		fprintf(['Loading mask volumes for subject number ' num2str(subjectNumber) ' of ' num2str(nSubjects) ': ' subject]);

		for maskNumber = 1:nMasks % For each mask...
	
			% Figure out which mask we're looking at
			mask = userOptions.maskNames{maskNumber};

			% Load the mask
			readPath = replaceWildcards(userOptions.maskPath, '[[subjectName]]', subject, '[[maskName]]', mask);
			maskMatrix = spm_read_vols(spm_vol(readPath));
			
			% Convert any NaNs to 0s (sometimes NaNs remain after masks are reverse-normalised
			maskMatrix(isnan(maskMatrix)) = 0;

			% Check to make sure it's not empty
			if numel(find(maskMatrix)) == 0, error(['The mask "' mask '" for subject "' thisSubject '" is empty!  Check your mask definitions.']); end%if
			
			% Store the loaded mask into a struct
			binaryMasks_nS.(subject).(mask) = maskMatrix;
			clear maskMatrix;
			
			fprintf('.');

		end%for:maskNumber
		
		fprintf(':\n');
		
	end%for:subjectNumber

	%% Save relevant info

	timeStamp = datestr(now);

	fprintf(['Saving masks to ' fullfile(userOptions.rootPath, 'ImageData', MasksFilename) '\n']);
	gotoDir(userOptions.rootPath, 'ImageData');
	save(MasksFilename, 'binaryMasks_nS','-v7.3');
	
	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');
	
else
	fprintf(['Loading previously saved masks from ' fullfile(userOptions.rootPath, 'ImageData', MasksFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'ImageData', MasksFilename));
end%if:overwriteFlag

if nargout == 1
	varargout{1} = binaryMasks_nS;
elseif nargout > 0
	error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started

end%function
