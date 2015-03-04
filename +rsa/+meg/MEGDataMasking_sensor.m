% MEGDataMasking_sensor is a function which takes a full set of MEG brain scans
% in sensor space and masks them according the the masks specified by the user.
%
% [maskedSensors =] MEGDataMasking_sensor(sensorImages, ...
%                                         maskSpec, ...
%                                         betaCorrespondence, ...
%                                         userOptions ...
%                                        )
%
%        sensorImages --- The unmasked sensor topographies.
%               sensorImages.(subjectName).(Grad1|Grad2|Mag|EEG) is a
%               [nChannels nTimepoints nConditions nSessions]-sized matrix.
%
%        maskSpec --- The specifications of the masking.
%               maskSpec.MEGSensorSites
%                       A vector containing the numbers of the MEG sensors to
%                       include. Defaults to all of them.
%               maskSpec.EEGSensorSites
%                       A vector containing the numbers of the EEG sensors to
%                       include. Leave blank to disregard EEG (untested!).
%                       Defaults to all of them.
%               maskSpec.baselineWindow
%                       A 2-sized vector containing the beginning and end of the
%                       time window of baseline, against which each sensor is
%                       normalised. The numbers refer to timepoints. Defaults to
%                       the first 50 timepoints ("[1 50]") if there are 50
%                       timepoints, otherwise the first n-1.
%               maskSpec.timeWindow
%                       A 2-sized vector containing the beginning and end of the
%                       time window of interest. The numbers refer to
%                       timepoints.Defaults to the whole epoch (not including
%                       baseline period). Woe betide those who let this overlap
%                       with the baselineWindow!
%               maskSpec.MEGSensors
%                       maskSpec.MEGSensors.Gradiometers
%                               Numerical. Can take the following values:
%                                       0  ---  Don't use gradiometers.
%                                       1  ---  Use gradiometers separately.
%                                       2  ---  Use gradiometers RMS.
%                               Defaults to 1.
%                       maskSpec.MEGSensors.Magnetometers
%                               Boolean. Use magnetometers? Defaults to false.
%                       maskSpec.MEGSensors.EEG
%                               Boolean. Use EEG? Defaults to false.
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
%
% The following files are saved by this function:
%        userOptions.rootPath/ImageData/
%                userOptions.analysisName_MaskedSensors.mat
%                        Contains the raw timecourses in a structure such that
%                        maskedSensors.(subjectName) is a [nChannels
%                        nConditions nSessions]-sized matrix (median across time
%                        window).
%        userOptions.rootPath/Details/
%                userOptions.analysisName_MEGDataMasking_sensor_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%  
% Cai Wingfield 12-2009, 6-2010, 7-2010

function [varargout] = MEGDataMasking_sensor(sensorImages, userOptions)

returnHere = pwd; % We'll return to the pwd when the function has finished

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('MEGDataMasking_sensor:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('MEGDataMasking_sensor:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(sensorImages));

maskSpec = userOptions.maskSpec;
betaCorrespondence = userOptions.betaCorrespondence;
%% %%
%% Set defaults for maskSpec
%% %%

% Not all sensors may be used, so only use the ones in the first subject (assume they're all the same)
availableSensorTypes = fieldnames(sensorImages.(userOptions.subjectNames{1}));

% Go through each sensor type and figure out whether we're using MEG and/or EEG
useMEG = false;
useEEG = false;
exampleSensorTypeUsed = '';
for sensorType = 1:numel(availableSensorTypes)
	switch availableSensorTypes{sensorType}
		case {'Grad1', 'Grad2', 'Mag'}
			useMEG = true;
            exampleSensorTypeUsed = availableSensorTypes{sensorType};
		case 'EEG'
			useEEG = true;
            exampleSensorTypeUsed = availableSensorTypes{sensorType};
	end%switch
end%for
if strcmpi(exampleSensorTypeUsed, ''), error('MEGDataMasking_sensor:NoSensors', '(!) No sensors in data!'); end%if

% Get numbers of channels
if useMEG
	nMEGChannels = size(sensorImages.(userOptions.subjectNames{1}).Grad1, 1);
else
	nMEGChannels = 0;
end%if
if useEEG
	nEEGChannels = size(sensorImages.(userOptions.subjectNames{1}).EEG, 1);
else
	nEEGChannels = 0;
end%if

% Set up maskSpec
% maskSpec = setIfUnset(maskSpec, 'MEGSensorSites', 1:nMEGChannels);
% maskSpec = setIfUnset(maskSpec, 'EEGSensorSites', 1:nEEGChannels);
maskSpec = setIfUnset(maskSpec, 'baselineWindow', [1 min([50, size(sensorImages.(userOptions.subjectNames{1}).(exampleSensorTypeUsed)) - 1])]);
% maskSpec = setIfUnset(maskSpec, 'timeWindow', [maskSpec.baselineWindow(2)+1 size(sensorImages.(userOptions.subjectNames{1}).(exampleSensorTypeUsed), 2)]);
% maskSpec = setIfUnset(maskSpec, 'MEGSensors', struct());
% maskSpec.MEGSensors = setIfUnset(maskSpec.MEGSensors, 'Gradiometers', useMEG);
% maskSpec.MEGSensors = setIfUnset(maskSpec.MEGSensors, 'Magnetometers', useMEG);
% maskSpec.MEGSensors = setIfUnset(maskSpec.MEGSensors, 'EEG', useEEG);
% maskSpec = setIfUnset(maskSpec, 'patternType', 'Spatial');

%% %%
%% maskSpec done!
%% %%

% The filenames contain the analysisName as specified in the user options file
MaskedSensorsFilename = [userOptions.analysisName, '_MaskedSensors.mat'];
DetailsFilename = [userOptions.analysisName, '_MEGDataMasking_Details.mat'];

promptOptions.functionCaller = 'MEGDataMasking_sensor';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'ImageData', MaskedSensorsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
	
	nSubjects = numel(userOptions.subjectNames);
	
	betas = betaCorrespondence;
	
	nConditions = size(betas, 2);
	nSessions = size(betas,1);
	
	% Make the intervals in maskSpec into actual time-index lists
	baselineWindow = maskSpec.baselineWindow(1):maskSpec.baselineWindow(2);
	timeWindow = userOptions.maskSpec.toDataPoints(1):userOptions.maskSpec.toDataPoints(2);
    
    thisMask = maskSpec.maskName{1};

	for subject = 1:nSubjects % For each subject

		% Figure out the subject's name
		thisSubject = userOptions.subjectNames{subject};

		% Get the subjects data
		subjectSensorData = sensorImages.(thisSubject);

		% Which sensors?
		allSensorTypes = fieldnames(subjectSensorData);
		
		% Remove unwanted sensors:
		for sensorType = 1:numel(allSensorTypes) % For each sensor type
			thisSensorType = allSensorTypes{sensorType};
			switch thisSensorType
				% Remove them if they're not wanted (as specified in maskSpec)
				case 'Grad1'
					if ~maskSpec.MEGSensors.Gradiometers
						subjectSensorData = rmfield(subjectSensorData, 'Grad1');
					end%if
				case 'Grad2'
					if ~maskSpec.MEGSensors.Gradiometers
						subjectSensorData = rmfield(subjectSensorData, 'Grad2');
					end%if
				case 'Mag'
					if ~maskSpec.MEGSensors.Magnetometers
						subjectSensorData = rmfield(subjectSensorData, 'Mag');
					end%if
				case 'EEG'
					if ~maskSpec.MEGSensors.EEG
						subjectSensorData = rmfield(subjectSensorData, 'EEG');
					end%if
			end%switch
		end%for:sensorType
		
		if maskSpec.MEGSensors.Gradiometers
			subjectSensorData.RMS = zeros(size(subjectSensorData.Grad1));
		end%if:RMS
		
        for session = 1:nSessions % For each session...
            for condition = 1:nConditions % and each condition...
                
                % RMS Grads if desired:
                if maskSpec.MEGSensors.Gradiometers
                    % If both grad channels are there
                    if isfield(subjectSensorData, 'Grad1') && isfield(subjectSensorData, 'Grad2')
                        if ~isempty(subjectSensorData.Grad1) && ~isempty(subjectSensorData.Grad2)
                            for sensor = 1:size(subjectSensorData.Grad1, 1)
                                subjectSensorData.RMS(sensor, :, condition, session) = sqrt((subjectSensorData.Grad1(sensor, :, condition, session).^2 + subjectSensorData.Grad2(sensor, :, condition, session).^2) ./ 2);
                            end%for:sensor
                        else
                            warning('One of the gradiometers does not have data. RMS not calculated.')
                            
                        end
                    else
                        warning('Not both gradiometer channels are present, RMS not calculated.');
                    end%if
                end%if:RMS
                
            end%for:condition
        end%for:session
		
		% RMS Grads if desired:
		if maskSpec.MEGSensors.Gradiometers == 2
			% Remove old non-RMS'd grad channels
			subjectSensorData = rmfield(subjectSensorData, 'Grad1');
			subjectSensorData = rmfield(subjectSensorData, 'Grad2');
		end%if
		
		% Remaining sensor types
		remainingSensorTypes = fieldnames(subjectSensorData);

		for session = 1:nSessions % For each session...
			for condition = 1:nConditions % and each condition...
				
				overallSensorNumber = 1;
			
				% Now do the masking
				for sensorType = 1:numel(remainingSensorTypes)

					% Which sensor?
					thisSensorType = remainingSensorTypes{sensorType};

					% Before the masking, Get the data for the baseline (for normalisation)
					baselineMatrix = zeros(size(subjectSensorData.(thisSensorType), 1), numel(baselineWindow));
					baselineMatrix(:,:) = squeeze(subjectSensorData.(thisSensorType)(:, baselineWindow, condition, session)); % [sensor, time]
					
					% Work out (chanel-)mean (time-)standard deviations for normalisation factors and normalise
					baselineStd = std(baselineMatrix, 0, 2); % [sensor 1]
					baselineMStds = mean(baselineStd, 1); % [1 1]
					baselineMStds = squeeze(baselineMStds); % [1]
					normalisedSubjectSensorData = subjectSensorData.(thisSensorType)(:, :, condition, session) / baselineMStds;
					
					% Mask out the unwanted sensors and mask in the time window of interest
					switch thisSensorType
						case {'Grad1', 'Grad2', 'Mag', 'RMS'}
							maskedNormalisedSubjectSensorData = normalisedSubjectSensorData(maskSpec.MEGSensorSites, timeWindow);
						case 'EEG'
							maskedNormalisedSubjectSensorData = normalisedSubjectSensorData(maskSpec.EEGSensorSites, timeWindow);
					end%switch:thisSensorType
					
					% Forget which sensor type is which
					for sensor = 1:size(maskedNormalisedSubjectSensorData, 1)
						agnosticSensorData(overallSensorNumber, :, condition, session) = maskedNormalisedSubjectSensorData(sensor, :);
						overallSensorNumber = overallSensorNumber + 1;
					end%for:sensor
				end%for:sensorType
			end%for:condition
		end%for:session
	
		% Transform the sensor data into the appropriate type of 'pattern':
		switch lower(maskSpec.patternType)
		
			% For spatial patterns, median across the time window
			case 'spatial'
				subjectMatrix = zeros(size(agnosticSensorData, 1), size(agnosticSensorData, 3), size(agnosticSensorData, 4));
				subjectMatrix(:,:,:) = median(agnosticSensorData, 2);
			
			% For temporal patterns, average across the sensors.
			case 'temporal'
				
				subjectMatrix = zeros(size(agnosticSensorData, 2), size(agnosticSensorData, 3), size(agnosticSensorData, 4));
				subjectMatrix(:,:,:) = mean(agnosticSensorData, 1);
				
			% For spatiotemporal patterns, concatenate the sensor timecourses
			case 'spatiotemporal'
				
				subjectMatrix = zeros(1, size(agnosticSensorData, 3), size(agnosticSensorData, 4));
				for sensor = 1:size(agnosticSensorData, 1)
					if sensor == 1
						subjectMatrix(:,:,:) = agnosticSensorData(sensor, :, :, :);
					else
						subjectMatrix = [subjectMatrix; agnosticSensorData(sensor, :, :, :)];
					end%if
				end%for:sensor
				
		end%switch:patternType

		% For each subject, record the vectorised brain scan in a subject-name-indexed structure
		maskedSensors.(thisMask).(thisSubject) = subjectMatrix;
		
		clear subjectSensorData normalisedSubjectSensorData agnosticSensorData;

	end%for:subject

	%% Save relevant info

	timeStamp = datestr(now);

	fprintf(['Saving masked data to ' fullfile(userOptions.rootPath, 'ImageData', MaskedSensorsFilename) '\n']);
	save(fullfile(userOptions.rootPath, 'ImageData', MaskedSensorsFilename), 'maskedSensors');
	
	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	save(fullfile(userOptions.rootPath, 'Details', DetailsFilename), 'timeStamp', 'userOptions', 'maskSpec');
	
else
	fprintf(['Loading previously saved masked topographies from ' fullfile(userOptions.rootPath, 'ImageData', MaskedSensorsFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'ImageData', MaskedSensorsFilename));
end%if

if nargout == 1
	varargout{1} = maskedSensors;
elseif nargout > 0
	error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere); % Go back