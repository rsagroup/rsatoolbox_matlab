% fMRISearchlight is a function which takes some full brain volumes of data,
% some binary masks and some models and perfoms a searchlight in the data within
% each mask, matching to each of the models.  Saved are native-space r-maps for
% each model.
%
% [rTopographies, pTopographies, searchlightRDMs[, nTopographies]
%                                   =] MEGSearchlight_sensor(
%                                                            sensorImages,
%                                                            Models,
%                                                            userOptions
%                                                           )
%
%        sensorImages --- The unmasked sensor topographies.
%               sensorImages.(subjectName).(Grad1|Grad2|Mag|EEG) is a
%               [nChannels nTimepoints nConditions nSessions]-sized matrix.
%
%        Models --- A stack of model RDMs in a structure.
%               models is a [1 nModels] structure with fields:
%                       RDM
%                       name
%                       color
%
%        maskSpec --- The specifications of the masking.
%               maskSpec.timeWindow
%                       A 2-sized vector containing the beginning and end of the
%                       time window of interest. The numbers refer to
%                       timepoints. Defaults to the whole epoch.
%               maskSpec.MEGSensors
%                       maskSpec.MEGSensors.Gradiometers
%                               Numerical. Can take the following values:
%                                       0  ---  Don't use gradiometers.
%                                       1  ---  Use gradiometers separately.
%                                       2  ---  Use gradiometers RMS.
%                               Defaults to 1.
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
%                        names. Defaults to the fieldnames in sensorImages.
%                userOptions.conditionLabels
%                        A cell array containing the names of the conditions in
%                        this experiment.
%                userOptions.distance
%                        A string indicating the distance measure with which to
%                        calculate the RDMs. Defaults to 'Correlation'.
%                userOptions.distanceMeasure
%                        A string descriptive of the distance measure to be used
%                        to compare two RDMs. Defaults to 'Spearman'.
%                userOptions.temporalSearchlightResolution
%                        An integer. The number of timepoints which the temporal
%                        sliding window will move by for each iteration of the
%                        searchlight. Defaults to 1.
%                userOptions.temporalSearchlightWidth
%                        An integer. The width of the temporal sliding window,
%                        in timepoints. Defaults to 1.
%                userOptions.sensorSearchlightRadius
%                        An integer. The number of steps away from the centre of
%                        the searchlight that it extends. 0 means there's only
%                        one sensor. 1 means there's the centre, plus everything
%                        adjacent to it. 2 means everything 2-adjacent to the
%                        centre, etc. Defaults to 1.
%
% The following files are saved by this function:
%        userOptions.rootPath/Maps/
%                userOptions.analysisName_MEGSensorSearchlight_Maps.mat
%                        Contains the searchlight statistical maps in struct
%                        containing rTopogrpahies, pTopographies and
%                        nTopogrhies.
%                Various .fif files containing the searchlight results.
%                        >>> A NOTE ABOUT THESE FILES:
%                                Due to a bug/feature in fiff_write_evoked, it
%                                is necessary to have the first time point of a
%                                fiff file to be labled with a negative index.
%                                Because this won't work in every scenario,
%                                Instead here the times are ARBITRARY, so that 0
%                                is roughly in the middle.  A time index should
%                                only be considered as relitive to the first one
%                                and referenced against the user preference for
%                                this. Sorry!
%        userOptions.rootPath/RDMs/
%                userOptions.analysisName_MEGSensorSearchlight_RDMs.mat
%                        Contains the RDMs for each searchlight so that
%                        searchlightRDMs.(subject)(:, :, x, y, z) is the RDM.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_MEGSensorSearchlight_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
% Cai Wingfield 3-2010 7-2010
% Edited by IZ 07/13 

function [varargout] = MEGSearchlight_sensor(sensorImages, Models, userOptions)

	returnHere = pwd; % We'll come back here later

	%% Set defaults and check options struct
	if ~isfield(userOptions, 'analysisName'), error('MEGSearchlight_sensor:NoAnalysisName', 'analysisName must be set. See help'); end%if
	if ~isfield(userOptions, 'rootPath'), error('MEGSearchlight_sensor:NoRootPath', 'rootPath must be set. See help'); end%if
	userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(sensorImages));
	userOptions = setIfUnset(userOptions, 'temporalSearchlightResolution', 1);
	userOptions = setIfUnset(userOptions, 'temporalSearchlightWidth', 1);
	userOptions = setIfUnset(userOptions, 'sensorSearchlightRadius', 1);
	userOptions = setIfUnset(userOptions, 'distanceMeasure', 'Spearman');
	userOptions = setIfUnset(userOptions, 'distance', 'Correlation');
	
    betaCorrespondence = userOptions.betaCorrespondence;
	
	% The analysisName will be used to label the files which are eventually saved.
	MapsFilename = [userOptions.analysisName, '_MEGSensorSearchlight_Maps.mat'];
	RDMsFilename = [userOptions.analysisName, '_MEGSensorSearchlight_RDMs.mat'];
	DetailsFilename = [userOptions.analysisName, '_MEGSensorSearchlight_Details.mat'];

	promptOptions.functionCaller = 'MEGSearchlight_sensor';
	promptOptions.defaultResponse = 'S';
	promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', MapsFilename);
	promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

	overwriteFlag = overwritePrompt(userOptions, promptOptions);

	if overwriteFlag
		
		fprintf('Shining RSA searchlights...\n');
		
		% Data
		nSubjects = numel(userOptions.subjectNames);
		
		% Remove EEG if it's there
		if isfield(sensorImages.(userOptions.subjectNames{1}), 'EEG')
			fprintf('        Sensor searchlight not implemented for EEG, these data will be ignored.\n');
			for subject = 1:nSubjects
				sensorImages.(userOptions.subjectNames{subject}) = rmfield(sensorImages.(userOptions.subjectNames{subject}), 'EEG');
			end%for:subject
		end%if:EEG
        
        if isfield(sensorImages.(userOptions.subjectNames{1}), 'Mag')
            fprintf('        Sensor searchlight not implemented for Mag, these data will be ignored.\n');
            for subject = 1:nSubjects
                sensorImages.(userOptions.subjectNames{subject}) = rmfield(sensorImages.(userOptions.subjectNames{subject}), 'Mag');
            end%for:subject
        end%if:EEG
		
		% Do RMS if necessary
        if isfield(sensorImages.(userOptions.subjectNames{1}), 'Grad1') && isfield(sensorImages.(userOptions.subjectNames{1}), 'Grad2')
            for subject = 1:nSubjects
                
                % RMS
                sensorImages.(userOptions.subjectNames{subject}).RMS = sqrt((sensorImages.(userOptions.subjectNames{subject}).Grad1 .^2 + sensorImages.(userOptions.subjectNames{subject}).Grad2 .^2) ./ 2);
                
                % Remove Grads
                sensorImages.(userOptions.subjectNames{subject}) = rmfield(sensorImages.(userOptions.subjectNames{subject}), 'Grad1');
                sensorImages.(userOptions.subjectNames{subject}) = rmfield(sensorImages.(userOptions.subjectNames{subject}), 'Grad2');
                
            end%for:subject
        else
            warning('Not both gradiometer channels are present, RMS not calculated.');
        end%if:Grads
        
		
		% Which sensors do we have left?
		availableSensorTypes = fieldnames(sensorImages.(userOptions.subjectNames{1}));
        if numel(availableSensorTypes) == 0
            error('No gradiometers found! Please check your data.')
        end
		
		nSensorSites = size(sensorImages.(userOptions.subjectNames{1}).(availableSensorTypes{1}), 1);
		
		searchlightOptions.monitor = false;
		searchlightOptions.fisher = true;

		for subjectNumber = 1:nSubjects % and for each subject...
		
			tic;%1

			fprintf(['\tSearching in the sensor topographies of subject ' num2str(subjectNumber) ' of ' num2str(nSubjects)]);

			% Figure out which subject this is
			subject = userOptions.subjectNames{subjectNumber};
			
			tempBetas = betaCorrespondence;

			searchlightOptions.nSessions = size(tempBetas, 1);
			searchlightOptions.nConditions = max(size(userOptions.conditionLabels));
			searchlightOptions.sensorSearchlightPatterns = userOptions.searchlightPatterns;
			
			% Full brain data volume to perform searchlight on
			singleSubjectTopography = sensorImages.(subject);
			
			for sensorType = 1:numel(availableSensorTypes)
				
				% Which sensor?
				thisSensorType = availableSensorTypes{sensorType};
				
				%%## NORMALISE?!
				
				% Look only in relevant temporal limits
				singleSubjectTopography.(thisSensorType) = singleSubjectTopography.(thisSensorType)(:, userOptions.toDataPoints(1):userOptions.toDataPoints(2), :, :);
				
			end%for:sensorType

			% Do the searchlight!
			[rs, ps, ns, searchlightRDMs.(subject)] = searchlightMapping_MEG_sensor(singleSubjectTopography, Models, userOptions, searchlightOptions); % ps are from linear correlation p-values, and so aren't too useful here.

			nTopographies.(subject) = ns; % How many sensors contributed to the searchlight centred at each point. (Those with n==1 are excluded because the results aren't multivariate.)

			% Load metadata struct for this subject
			readFile = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', tempBetas(1,1).identifier, '[[subjectName]]', subject); %1,1 should be fine
			[ignore, subjectMetadataStruct] = evalc('fiff_read_evoked(readFile)');

            modelName = spacesToUnderscores(Models(userOptions.modelNumber).name);
            
            % Store results in indexed volumes
            rTopographies.(subject) = rs(:,:); % r-values for correlation with each model
            pTopographies.(subject) = ps(:,:); % p-values for correlation with each model
            
            % Parameters
            
            epochLength = size(singleSubjectTopography.(availableSensorTypes{1}), 2); % in time points
            nTimePoints = size(rs, 2);
            
            % Triplicate each sensor for writing to a .fif file (just in
            % case).
            fifRs = nan([3*nSensorSites, nTimePoints]);
            for sensorSite = 1:nSensorSites
                fifRs(sensorSite * 3 - 2, :) = squeeze(rs(sensorSite,:));
                fifRs(sensorSite * 3 - 1, :) = squeeze(rs(sensorSite,:));
                fifRs(sensorSite * 3, :) = squeeze(rs(sensorSite,:));
            end%for:sensorSites
            
            % sfreq is in samples per second
            % the original sfreq will be nTimePoints_orig / timeCourseLength_seconds
            % the new one will be nTimePoints / timeWindowLength_seconds
            sfreq = floor(subjectMetadataStruct.info.sfreq / userOptions.temporalSearchlightResolution); 
            subjectMetadataStruct.info.sfreq = sfreq; % sample frequency in samples per second?
            
            first = ceil(-nTimePoints / 2); % first time point (mne => this must be negative)
            last = floor(nTimePoints / 2) - 1; % last time point (minus 1... because I have to :( )
            times = linspace(first,last,nTimePoints)/sfreq;%(userOptions.temporalSearchlightLimits(1):userOptions.temporalSearchlightResolution:epochLength-userOptions.temporalSearchlightWidth) / 1000; % time in seconds
            
            % The unused channels, including EEG, are set to 0
            epochs = subjectMetadataStruct.evoked.epochs(:,1:size(fifRs,2));
            epochs(1:size(fifRs,1),:) = fifRs;
            
            % Put this all in a mne metadata struct
            subjectMetadataStruct.evoked.first = first;
            subjectMetadataStruct.evoked.last = last;
            subjectMetadataStruct.evoked.times = times;
            subjectMetadataStruct.evoked.epochs = epochs;
            
            % Write the fiff file
            outputFilename = [userOptions.analysisName '_rTopography_' modelName '_' subject '.fif'];
            if ~exist(fullfile(userOptions.rootPath, 'Maps'), 'dir')
                mkdir(fullfile(userOptions.rootPath, 'Maps'));
            end
            fiff_write_evoked(fullfile(userOptions.rootPath, 'Maps',outputFilename), subjectMetadataStruct);
            
            % Convert this fiff file to an SPM file UNFINISHED!
            %fifFileName = outputFilename;
            %fif2spm
            
            
            clear rs ps ns;

			t = toc;%1
			fprintf([': [' num2str(ceil(t)) 's]\n']);
			
		end%for:subjectNumber

		%% Save relevant info

		timeStamp = datestr(now);

		fprintf(['Saving searchlight maps to ' fullfile(userOptions.rootPath, 'Maps', MapsFilename) '\n']);
		gotoDir(userOptions.rootPath, 'Maps');
		save(MapsFilename, 'rTopographies', 'pTopographies', 'nTopographies');
		
		fprintf(['Saving RDMs to ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '\n']);
		gotoDir(userOptions.rootPath, 'RDMs');
		save(RDMsFilename, 'searchlightRDMs');
		
		fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
		gotoDir(userOptions.rootPath, 'Details');
		save(DetailsFilename, 'timeStamp', 'userOptions');
		
	else
		fprintf(['Loading previously saved maps from ' fullfile(userOptions.rootPath, 'Maps', MapsFilename) '...\n']);
		load(fullfile(userOptions.rootPath, 'Maps', MapsFilename));
	end%if

	if nargout == 3
		varargout{1} = rTopographies;
		varargout{2} = pTopographies;
		varargout{3} = searchlightRDMs;
	elseif nargout == 4
		varargout{1} = rTopographies;
		varargout{2} = pTopographies;
		varargout{3} = searchlightRDMs;
		varargout{4} = nTopographies;
	elseif nargout > 0
		error('0, 3 or 4 arguments out, please.');
	end%if:nargout

	cd(returnHere); % And go back to where you started

end%function
