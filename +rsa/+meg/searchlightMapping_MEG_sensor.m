%% %%
%% Subfunctions
%% %%

function [smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_MEG_sensor(singleSubjectTopography, models, userOptions, localOptions)

	%  [smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_MEG_sensor(singleSubjectTopography, models, mask, userOptions, localOptions)
	%
	% CW 5-2010 Edited IZ 07/13 

	%% Get parameters

	% Prepare models
    nModels = size(models,2);
    if nModels>1
        modelRDMs_utv = squeeze(unwrapRDMs(vectorizeRDMs(models)))'; % note by IZ: if Models has 1 model, this is a column vector, else models are put into rows
    else
        modelRDMs_utv = squeeze(unwrapRDMs(vectorizeRDMs(models)));
    end;

	availableSensorTypes = fieldnames(singleSubjectTopography);

	nConditions = size(singleSubjectTopography.(availableSensorTypes{1}), 3);
	nSessions = size(singleSubjectTopography.(availableSensorTypes{1}), 4);

	% Get the sensor layout map
	%adjacencyMatrixLocation = fullfile(userOptions.toolboxRoot, 'Engine', 'Cai', 'neuromagSensors.csv');
	adjacencyMatrix = importSensorLayout('neuromagSensors.csv');

	% How long is the stimulus (in data points)?
	epochLength = numel(singleSubjectTopography.(availableSensorTypes{1})(1,:,1,1)); % Just need to look at one sensor type and one condition/session/subject

	% Number of sensors
	nSensorSites = 102;
	nSensorTypes = numel(fieldnames(singleSubjectTopography));

	% Number of time points to loop with
	nTimePoints = floor((epochLength - userOptions.temporalSearchlightWidth) / userOptions.temporalSearchlightResolution);

	%% similarity-graph-map the volume with the searchlight

	% Preallocate looped matrices for speed
	smm_ps = nan([nSensorSites, nTimePoints, nModels]);
	smm_rs = nan([nSensorSites, nTimePoints, nModels]);
	n = nan(nSensorSites, 1);
	searchlightRDMs = nan([nConditions, nConditions, nSensorSites, nTimePoints]);

	for sensorSite = 1:nSensorSites

		% Determine which sensor sites are within the radius of the currently-picked sensor site
		sensorSitesWithinRadius = getNearbySensors(sensorSite, adjacencyMatrix, userOptions.sensorSearchlightRadius);
		
		% How many?
		nSensorSitesWithinRadius = numel(sensorSitesWithinRadius);
		n(sensorSite) = nSensorSitesWithinRadius;

		for t = 1:nTimePoints

			% Work out the current time window
			currentTimeStart = 1 + (t - 1) * userOptions.temporalSearchlightResolution;
			currentTimeWindow = (currentTimeStart : currentTimeStart + userOptions.temporalSearchlightWidth - 1);

			% Restrict data to current searchlight
			for sensorType = 1:numel(availableSensorTypes)
				thisSensorType = availableSensorTypes{sensorType};
				if sensorType == 1
					currentData = singleSubjectTopography.(thisSensorType)(sensorSitesWithinRadius, currentTimeWindow, :, :); % (sensors, time, condition, session)
				else
					currentData = cat(1, currentData, singleSubjectTopography.(thisSensorType)(sensorSitesWithinRadius, currentTimeWindow, :, :));
				end%if:firstSensorType
			end%for:sensorType

			% Median over the time window
			switch lower(localOptions.sensorSearchlightPatterns)
				case 'spatial'
					% Spatial patterns: median over time window
					currentData = median(currentData, 2);
                    currentData = squeeze(currentData);
				case 'temporal'
					% Temporal patterns: mean over sensors within searchlight
					currentData = mean(currentData, 1);
                    currentData = squeeze(currentData);
				case 'spatiotemporal'
					currentData = reshape(currentData, [], size(currentData, 3), size(currentData, 4)); % (dataPoints, conditions, sessions)

			end%switch:localOptions.sensorSearchlightPatterns

			searchlightRDM = zeros(nConditions, nConditions);
			
			for session = 1:localOptions.nSessions

				searchlightRDM = searchlightRDM + squareform(pdist(squeeze(currentData(:,:,:,session))', userOptions.distance));
			
			end%for:sessions
			
			searchlightRDM = searchlightRDM ./ localOptions.nSessions;
		
			searchlightRDM = vectorizeRDM(searchlightRDM);
			
			% Locally store the full brain's worth of indexed RDMs.
			searchlightRDMs(:,:,sensorSite, t) = squareform(searchlightRDM);
		
			try
				[rs, ps] = corr(searchlightRDM', modelRDMs_utv', 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
			catch
				[rs, ps] = corr(searchlightRDM', modelRDMs_utv, 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
			end%try
		
			if localOptions.fisher
				for i = 1:numel(rs)
					rs(i) = fisherTransform(rs(i));
				end%for:i
			end%if
			
			smm_ps(sensorSite, t, :) = ps;
			smm_rs(sensorSite, t, :) = rs;

		end%for:t
		
		% Indicate progress every once in a while...
		if mod(sensorSite, 10) == 0, fprintf('.'); end%if
		
	end%for:sensorSite

end%function