% MEGDataPreparation_sensor is a function designed to take MEG data and load it into
% the correct format for the rest of the toolbox to use.
%
% [sensorImages[, baselineLimit] =] MEGDataPreparation_sensor( ...
%                                                        betaCorrespondence, ...
%                                                        userOptions)
%
% Specifics:  All three sensor types are normalised by their mean baseline
% standard deviation (making them into comparable SNRs (thanks to Olaf!)).
% Medians are then taken across their time windows of interest to produce a
% single vector for each subject-condition-scan triple.
%
% It is based on Su Li's "LexPro_MEG" script of 12-2009
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
%                userOptions.analysisName_SensorImages.mat
%                        Contains the raw timecourses in a structure such that
%                        sensorImages.(subjectName).(Grad1|Grad2|Mag|EEG|All) is
%                        a [nChannels nTimepoints nConditions nSessions]-sized
%                        matrix.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_MEGDataPreparation_sensor_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
% Cai Wingfield 12-2009, 3-2010, 6-2010

function [varargout] = MEGDataPreparation_sensor(userOptions)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll return to the pwd when the function has finished

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('MEGDataPreparation_sensor:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('MEGDataPreparation_sensor:NoRootPath', 'rootPath must be set. See help'); end%if
if ~isfield(userOptions, 'betaPath'), error('MEGDataPreparation_sensor:NoBetaPath', 'betaPath must be set. See help'); end%if
if ~isfield(userOptions, 'subjectNames'), error('MEGDataPreparation_sensor:NoSubjectNames', 'betaPath must be set. See help'); end%if

% The filenames contain the analysisName as specified in the user options file
ImageDataFilename = [userOptions.analysisName, '_SensorImages.mat'];
DetailsFilename = [userOptions.analysisName, '_MEGDataPreparation_sensor_Details.mat'];
BaselineLimitFilename = [userOptions.analysisName, '_BaselineLimit.mat'];

promptOptions.functionCaller = 'MEGDataPreparation_sensor';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag

	%% Get Data

	betas = userOptions.betaCorrespondence;
	
	nConditions = size(betas, 2);
	nSessions = size(betas,1);
	nSubjects = numel(userOptions.subjectNames);
	
	% These don't ever change
	MEGChannels = 1:102; % All the MEG channels
	grad1Channels = 3 .* MEGChannels - 2;
	grad2Channels = grad1Channels + 1;
	magChannels = grad1Channels + 2;

	for subject = 1:nSubjects % For each subject

		% Figure out the subject's name
		thisSubject = userOptions.subjectNames{subject};

		fprintf(['Reading MEG time courses for subject number ' num2str(subject) ' of ' num2str(nSubjects) ': ' thisSubject ':']);

		for session = 1:nSessions % For each session...  UNUSED?!?!
			for condition = 1:nConditions % and each condition...

				% Then read the brain data (for this session, condition)
				readPath = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', betas(session, condition).identifier, '[[subjectName]]', thisSubject);
				[ignore, allMEGData] = evalc('fiff_read_evoked(readPath)'); % Using evalc supresses output!
				
				nChannels = numel(allMEGData.info.chs);
				[EEGChannelsCell{1:nChannels}] = deal(allMEGData.info.chs.kind);
				EEGChannels = cell2mat(EEGChannelsCell);
				clear nChannels EEGChannelsCell;
				EEGChannelsBin = EEGChannels == 2;

				allMEGDataMatrix = allMEGData.evoked.epochs;
				
				baselineLimit = double(-allMEGData.evoked.first); % I hope these are all the same!

				% Get the data for the time window of interest
				grad1DataMatrix = allMEGDataMatrix(grad1Channels, :);
				grad2DataMatrix = allMEGDataMatrix(grad2Channels, :);

				magDataMatrix = allMEGDataMatrix(magChannels, :);

				EEGDataMatrix = allMEGDataMatrix(EEGChannelsBin, :);

				subjectSensorData.Grad1(:, :, condition, session) = grad1DataMatrix; % (channel, time, condition, session)
				subjectSensorData.Grad2(:, :, condition, session) = grad2DataMatrix;
				subjectSensorData.Mag(:, :, condition, session) = magDataMatrix;
				subjectSensorData.EEG(:, :, condition, session) = EEGDataMatrix;
				
				fprintf('\b.:');

			end%for
			
			fprintf('\n');

		end%for

		% For each subject, record the vectorised brain scan in a subject-name-indexed structure
		sensorImages.(thisSubject) = subjectSensorData;
		clear subjectSensorData;

	end%for

	%% Save relevant info

	timeStamp = datestr(now);

	fprintf(['Saving image data to ' fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename) '\n']);
	gotoDir(userOptions.rootPath, 'ImageData');
	save(ImageDataFilename, 'sensorImages');
	save(BaselineLimitFilename, 'baselineLimit');
	
	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');
	
else
	fprintf(['Loading previously saved topographies from ' fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'ImageData', ImageDataFilename));
	fprintf(['Loading previously saved baseline limit from ' fullfile(userOptions.rootPath, 'ImageData', BaselineLimitFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'ImageData', BaselineLimitFilename));
end%if

if nargout == 1
	varargout{1} = sensorImages;
elseif nargout == 2
	varargout{1} = sensorImages;
	varargout{2} = baselineLimit;
elseif nargout > 0
	error('0, 1 or 2 arguments out, please.');
end%if:nargout

cd(returnHere); % Go back
