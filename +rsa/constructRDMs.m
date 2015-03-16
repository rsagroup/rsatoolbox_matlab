function [varargout] = constructRDMs(responsePatterns, betaCorrespondence, userOptions)
%
% [RDMs =] constructRDMs(responsePatterns, betaCorrespondence, userOptions)
%
% constructRDMs is a function which takes a matrix of imaging data and
% computes RDMs from it, one RDM for each ROI, subject, and scanning
% session. The RDMs are returned in a structure
% with names reflecting the ROI, subject and session. The returned RDMs are
% also saved. The RDM structure format ensures that the RDMs will be
% properly labeled in subsequent analyses.
%
% responsePatterns: The input response patterns. It contains the masked
% brain activity maps in a structure such that
% responsePatterns.(mask).(subject) is a [nMaskedVoxels nConditions
% nSessions] matrix, where mask is a binary brain map that defines the ROI.
% 
% betaCorrespondence: The array of beta filenames.
% betas(condition, session).identifier is a string which refers to the
% filename (not including the path) of the SPM beta image. Alternatively,
% this can be the string "SPM", in which case the SPM metadata will be used
% to infer this information, provided that userOptions.conditionLabels is
% set, and the condition labels are the same as those used in SPM. It must
% be noted that this module assumes that data is in such a format that the
% input structure has one field per subject, session and mask. Defining
% betaCorrespondence or choosing ‘SPM’ does not automatically compute the
% responsePatterns structure. At this point only the number of sessions
% and/or conditions are extracted from SPM or betaCorrespondence.
% 
% userOptions: The standard options structure containing the following
% fields:
% userOptions.analysisName: A string which is prepended to the  saved files.
% userOptions.rootPath: A string describing the root path  where files will
% be saved (inside created directories).
% userOptions.maskNames: A cell array containing strings identifying the
% mask names. It defaults to the field names of the first subject of
% responsePatterns. userOptions.subjectNames: A cell array containing
% strings identifying the subject names. Defaults to the field names in
% (the first mask of) responsePatterns.
% userOptions.distance: A string indicating the distance measure with which
% to calculate the RDMs. Defaults to “Correlation”, but can be set to any 
% Matlab distance measure.
% userOptions.RoIColor: A triple indicating the [R G B] value of the colour
% which should be used to indicated RoI RDMs on various diagrams. Defaults
% to black ([0 0 0]).
% The following files are saved by this function:
% userOptions.rootPath/RDMs/userOptions.analysisName_RDMs.mat Contains a
% structure of ROI RDMs which is of size [nMasks, nSubjects, nSessions] and
% with fields:, RDM, name, color. userOptions.rootPath/Details/
% userOptions.analysisName_constructRDMs_Details.mat is a file that
% contains the userOptions for this execution of the function and a
% timestamp.

% Cai Wingfield 11-2009, 12-2009, 3-2010, 6-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

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
if ~isfield(userOptions, 'analysisName'), error('constructRDMs:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('constructRDMs:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'maskNames', fieldnames(responsePatterns));
userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(responsePatterns.(userOptions.maskNames{1})));
userOptions = setIfUnset(userOptions, 'distance', 'Correlation');
userOptions = setIfUnset(userOptions, 'RoIColor', [0 0 0]);

% The analysisName will be used to label the files which are eventually saved.
RDMsFilename = [userOptions.analysisName, '_RDMs.mat'];
DetailsFilename = [userOptions.analysisName, '_constructRDMs_Details.mat'];

promptOptions.functionCaller = 'constructRDMs';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'RDMs', RDMsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
	
	%% Get Data

	if ischar(betaCorrespondence) && strcmpi(betaCorrespondence, 'SPM')
		betas = getDataFromSPM(userOptions);
	else
		betas = betaCorrespondence;
	end%if:SPM
	
	nSubjects = numel(userOptions.subjectNames);
	nMasks = numel(userOptions.maskNames);
	nSessions = size(betas, 1);
	nConditions = size(betas, 2);
	
	for mask = 1:nMasks % For each mask...

		thisMask = userOptions.maskNames{mask};

		for subject = 1:nSubjects % and for each subject...

			% Figure out which subject this is
			thisSubject = userOptions.subjectNames{subject};

			for session = 1:nSessions % and each session...

				% Get the brain scan vol
				thisActivityPattern = responsePatterns.(thisMask).(thisSubject);

				% Calculate the RDM
				localRDM = squareform( ...
					pdist( ...
						squeeze(thisActivityPattern(:, :, session))', userOptions.distance));

				% Store the RDM in a struct with the right names and things!
				RDMs(mask, subject, session).RDM = localRDM;
				RDMs(mask, subject, session).name = [deunderscore(thisMask) ' | ' thisSubject ' | Session: ' num2str(session)];
				RDMs(mask, subject, session).color = userOptions.RoIColor;

				clear localRDM;

			end%for:session
		end%for:subject
	end%for:mask

	%% Save relevant info

	timeStamp = datestr(now);

% 	fprintf(['Saving RDMs to ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '\n']);
    disp(['Saving RDMs to ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename)]);
	gotoDir(userOptions.rootPath, 'RDMs');
	save(RDMsFilename, 'RDMs');
	
% 	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
    disp(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename)]);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');
	
else
% 	fprintf(['Loading previously saved RDMs from ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '...\n']);
  	disp(['Loading previously saved RDMs from ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '...']);
    load(fullfile(userOptions.rootPath, 'RDMs', RDMsFilename));
end%if

if nargout == 1
	varargout{1} = RDMs;
elseif nargout > 0
	error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started
