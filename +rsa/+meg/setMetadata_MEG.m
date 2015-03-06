% userOptions = setMetadata_MEG(Models, userOptions)
%
% Models - an array of structs
% userOptions - struct of user preferences, modified and returned.
%
% This function sets any missing parameters in userOptions to the default
% values. If no default values can be set, ask the user to set it and show
% an error message. It also reads in the MEG data to find out some meta
% data related to the experiment.
%
% Li Su 3-2012
% Update: Isma Zulfiqar 9-2012 fixed time issues for searchlight all brain
% Update: IZ 11-2012 time support for ROI maks

function userOptions = setMetadata_MEG(Models, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if ~isfield(userOptions, 'analysisName'), error('projectOptions:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('projectOptions:NoRootPath', 'rootPath must be set. See help'); end%if
if ~isfield(userOptions, 'betaPath'), error('projectOptions:NoBetaPath', 'betaPath must be set. See help'); end%if
if ~isfield(userOptions, 'subjectNames'), error('projectOptions:NoSubjectNames', 'betaPath must be set. See help'); end%if

% loading one subject's STC/FIFFfile in order to find out its meta data.
tempBetas = betaCorrespondence();
userOptions.betaCorrespondence = tempBetas;

if userOptions.sensorLevelAnalysis
    
    % Despite warning, readPath *is* used, but it's used in an eval
    readPath = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', tempBetas(1,1).identifier, '[[subjectName]]', userOptions.subjectNames{1});
    
    %
    % Using evalc supresses output!
    [ignore, allMEGData] = evalc('fiff_read_evoked(readPath)');
    
    samplingRate = allMEGData.info.sfreq/1000; % ms
    
    tmin = double(allMEGData.evoked.first); tmax = double(allMEGData.evoked.last);
    
    if userOptions.searchlight
        time = userOptions.temporalSearchlightLimits;
        
        if time(1) < tmin
            disp(['Lower time limit is out of bounds. Using minimum value from fiff file instead. tmin = ' num2str(tmin) ' ms']);
            userOptions.temporalSearchlightLimits(1) = tmin;
        end
        if time(2) > tmax
            disp(['Upper time limit is out of bounds. Using maximum value from fiff file instead. tmax = ' num2str(tmax) ' ms']);
            userOptions.temporalSearchlightLimits(2) = tmax;
        end
        
        userOptions.toDataPoints = [1+((userOptions.temporalSearchlightLimits(1) - tmin)/samplingRate) 1+((userOptions.temporalSearchlightLimits(2)-tmin)/samplingRate)];
        
    else
        time = userOptions.maskSpec.timeWindow;
        
        if time(1) < tmin
            disp(['Lower time limit is out of bounds. Using minimum value from fiff file instead. tmin = ' num2str(tmin) ' ms']);
            userOptions.maskSpec.timeWindow(1) = tmin;
        end
        if time(2) > tmax
            disp(['Upper time limit is out of bounds. Using maximum value from fiff file instead. tmax = ' num2str(tmax) ' ms']);
            userOptions.maskSpec.timeWindow(2) = tmax;
        end
        
        userOptions.maskSpec.toDataPoints = [1+((userOptions.maskSpec.timeWindow(1) - tmin)/samplingRate) 1+((userOptions.maskSpec.timeWindow(2)-tmin)/samplingRate)];
        
    end
    
else % source level analysis
    
    % We can read just the left hemisphere, as we just want the metadata,
    % and we assume both hemispheres will be the same
    readPathL = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', tempBetas(1, 1).identifier, '[[subjectName]]', userOptions.subjectNames{1}, '[[LR]]', 'l');
    MEGDataStcL = mne_read_stc_file1(readPathL);
    MEGDataVolL = single(MEGDataStcL.data);
    
    userOptions.nSubjects = numel(userOptions.subjectNames);
    userOptions.monitor = false;
    userOptions.fisher = true;
    
    % TODO: This is bad
    userOptions.nVertices = userOptions.targetResolution; % Same for all subjects and both hemispheres, hopefully...
    
    % TODO: this is bad
    userOptions.nSessions = size(tempBetas, 1);
    modelNumber = userOptions.modelNumber;
    
    % Todo: this is bad
    userOptions.nConditions = size(squareRDM(Models(modelNumber).RDM), 1);
    
    % ============= setting time parameters for output file ================= %
    userOptions.STCmetaData.tmin = MEGDataStcL.tmin; % - ...
    % (MEGDataStcL.tmin - (userOptions.temporalSearchlightLimits(1) /1000)); % in seconds
    userOptions.STCmetaData.vertices = 1:userOptions.nVertices;
    userOptions.STCmetaData.tstep = MEGDataStcL.tstep;
    
    % time values are converted to their data point equivalents here.
    % this depends on the size of recordings input
    
    % ================= converting from mseconds to datapoints ============== %
    
    [vertices, totalDataPoints] = size(MEGDataStcL.data);
    totalTimeInMs = totalDataPoints * MEGDataStcL.tstep*1000; % calculated from time step and total data points
    timeAdjusted = totalTimeInMs - norm(MEGDataStcL.tmin * 1000); % last point of data from tmin and total time
    
    %============================= input checks ============================= %
    % ====== comparing search light resolution to the time step of data ===== %
    
    % time step
    if userOptions.searchlight % added 03/12 IZ
        if MEGDataStcL.tstep * 1000 * userOptions.temporalDownsampleRate > userOptions.temporalSearchlightResolution
            error('Error: The input time resolution of search light cannot be applied. The time resolution for data is lower.');
        else
            userOptions.STCmetaData.tstep = userOptions.temporalSearchlightResolution * ...
                userOptions.temporalDownsampleRate / 1000;  % in s
        end
    end
    
    % searchlight width in relation to downsample rate
    if userOptions.temporalSearchlightWidth < userOptions.STCmetaData.tstep * 1000
        disp('Warning: The searchlight width is less than data rate after downsampling');
        disp(['Setting temporal searchlight width equal to minimum timestep: ', num2str(userOptions.STCmetaData.tstep*1000), 'ms']);
        userOptions.temporalSearchlightWidth = userOptions.STCmetaData.tstep*1000;
    end
    
    %% MASKING TIME WINDOWS %%
    if userOptions.maskingFlag && not(userOptions.searchlight) % sliding time window and roi analysis
        
        nMasks = numel(userOptions.maskNames);
        for mask = 1:nMasks
            
            % Which mask is this?
            thisMask = dashToUnderscores(userOptions.maskNames{mask});
            
            time = cell2mat(userOptions.maskTimeWindows(mask));
            lowerTimeLimit = time(1);
            upperTimeLimit = time(2);
            
            % getting the starting data point by comparing starting points of searchlight and input data
            differenceInms = lowerTimeLimit - MEGDataStcL.tmin*1000;
            startingDataPoint = 1 + floor((differenceInms / (MEGDataStcL.tstep*1000)));
            
            % getting last data point for searchlight limits
            differenceInms = norm(upperTimeLimit - lowerTimeLimit);
            lastDataPoint = startingDataPoint + floor((differenceInms / (MEGDataStcL.tstep*1000)));
            
            % parameters for output file
            userOptions.STCmetaData.tmin = lowerTimeLimit /1000; % in seconds
            
            % checks onmask timing information
            if lowerTimeLimit < MEGDataStcL.tmin*1000
                disp(['Warning: The searchlight for mask: ', thisMask, ' is attempting to access time point not in data.'] );
                disp(strcat(' >> Using minimum time point value instead... ', int2str(MEGDataStcL.tmin*1000), ' ms'));
                startingDataPoint = 1;
                userOptions.STCmetaData.tmin = MEGDataStcL.tmin;
            end
            
            if userOptions.slidingTimeWindow % added by IZ 11-12
                if timeAdjusted-userOptions.temporalSearchlightWidth <= upperTimeLimit
                    disp(['Warning: The sliding window for mask: ', thisMask, ' over-runs the data sample. ']);
                    disp(strcat('>> Using maximum time point, with adjusted factor of source searchlight radius, from data instead... ', int2str(timeAdjusted-userOptions.sourceSearchlightRadius), ' ms'));
                    lastDataPoint = totalDataPoints - ...
                        ceil(userOptions.temporalSearchlightWidth/userOptions.temporalSearchlightResolution);
                end
            else
                if timeAdjusted <= upperTimeLimit
                    disp(['Warning: The searchlight for mask: ', thisMask, ' over-runs the data sample. ']);
                    disp(strcat('>> Using maximum time point from data instead... ', int2str(timeAdjusted-userOptions.sourceSearchlightRadius), ' ms'));
                    lastDataPoint = totalDataPoints;
                end
            end
            
            userOptions.maskTimetoDataPoints.(thisMask) = [startingDataPoint lastDataPoint];
        end
        
    elseif userOptions.searchlight % all brain searchlight
        
        userOptions.STCmetaData.tmin = userOptions.temporalSearchlightLimits(1)/1000;
        
        % getting the starting data point by comparing starting points of searchlight and input data
        differenceInms = userOptions.temporalSearchlightLimits(1) - MEGDataStcL.tmin*1000;
        startingDataPoint = 1 + floor((differenceInms / (MEGDataStcL.tstep*1000)));
        
        % getting last data point for searchlight limits
        differenceInms = norm(userOptions.temporalSearchlightLimits(2) - userOptions.temporalSearchlightLimits(1));
        lastDataPoint = startingDataPoint + floor((differenceInms / (MEGDataStcL.tstep*1000)));
        
        % lower temporal searchlight limit
        if userOptions.temporalSearchlightLimits(1) < MEGDataStcL.tmin*1000
            disp('Warning: The searchlight is attempting to access time point not in data.' );
            disp(strcat(' >> Using minimum time point value instead... ', int2str(MEGDataStcL.tmin*1000), ' ms'));
            startingDataPoint = 1;
            userOptions.STCmetaData.tmin = MEGDataStcL.tmin;
        end
        
        % upper temporal searchlight limit
        if timeAdjusted-userOptions.temporalSearchlightWidth <= userOptions.temporalSearchlightLimits(2)
            disp('Warning: The search light over-runs the data sample. ');
            disp(strcat('>> Using maximum time point, with adjusted factor of source searchlight radius, from data instead... ', int2str(timeAdjusted-userOptions.sourceSearchlightRadius), ' ms'));
            userOptions.temporalSearchlightLimits(2) = timeAdjusted;
            lastDataPoint = totalDataPoints - ...
                (userOptions.temporalSearchlightWidth/userOptions.temporalSearchlightResolution);
        end
        
        userOptions.dataPointsSearchlightLimits = [startingDataPoint lastDataPoint];
        
    end
    
    % ========================== new fields ================================ %
    userOptions.toDataPoints = totalDataPoints/(userOptions.temporalDownsampleRate*totalTimeInMs);
    
    % ====================================================================== %
    userOptions = setIfUnset(userOptions, 'minDist', 5);
    userOptions = setIfUnset(userOptions, 'maskingFlag', false);
    userOptions = setIfUnset(userOptions, 'primaryThreshold', 0.05);
    userOptions = setIfUnset(userOptions, 'ModelColor', [0 1 0]);
    userOptions = setIfUnset(userOptions, 'RoIColor', [0 0 1]);
    userOptions = setIfUnset(userOptions, 'temporalDownsampleRate', 1);
    userOptions = setIfUnset(userOptions, 'temporalSearchlightResolution', 1);
    userOptions = setIfUnset(userOptions, 'temporalSearchlightWidth', 1);
    userOptions = setIfUnset(userOptions, 'distanceMeasure', 'Spearman');
    userOptions = setIfUnset(userOptions, 'distance', 'Correlation');
    userOptions = setIfUnset(userOptions, 'sensorSearchlightRadius', 1);
    userOptions = setIfUnset(userOptions, 'sourceSearchlightRadius', userOptions.minDist);
    userOptions = setIfUnset(userOptions, 'temporalSearchlightLimits', [1 size(MEGDataVolL, 2)]);
    userOptions = setIfUnset(userOptions, 'searchlightPatterns', 'spatiotemporal');
    userOptions = setIfUnset(userOptions, 'significanceTestPermutations', 1000);
    userOptions = setIfUnset(userOptions, 'nResamplings', 1000);
    userOptions = setIfUnset(userOptions, 'groupStats', 'RFX');
    userOptions = setIfUnset(userOptions, 'targetResolution', 10242);
    userOptions = setIfUnset(userOptions, 'smoothingWidth', 10);
end
