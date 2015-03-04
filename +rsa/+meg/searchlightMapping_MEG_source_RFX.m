function [smm_rs, smm_ps] = searchlightMapping_MEG_source_RFX (singleMesh, Models, userOptions)

	%  [smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_MEG_source(singleMesh, Models, mask, userOptions, localOptions)
	%  based on Li Su's script
	% CW 5-2010, last updated by Li Su 3-2012

	%% Get parameters

	% Prepare Models
	modelRDMs_utv = squeeze(unwrapRDMs(vectorizeRDMs(Models)))';
	nModels=size(modelRDMs_utv,1);

	nVertices = userOptions.nVertices;
	nConditions = userOptions.nConditions;
	%nSessions = userOptions.nSessions;

	% How long is the stimulus (in time points)?
	epochLength = size(singleMesh, 2); % (vertex, time, condition, session)

	% Number of time points to loop with given the withd and time step of searchlight
	nTimePoints = floor((epochLength - userOptions.temporalSearchlightWidth) / userOptions.temporalSearchlightResolution);

	%% similarity-graph-map the volume with the searchlight

	% Preallocate looped matrices for speed
	smm_ps = nan([nVertices, nTimePoints, nModels]);
	smm_rs = nan([nVertices, nTimePoints, nModels]);
	%n = nan(nVertices);
	%searchlightRDMs = nan([nConditions, nConditions, nVertices, nTimePoints]);
    adjacencyMatrix = calculateMeshAdjacency(userOptions.nVertices, userOptions.sourceSearchlightRadius, userOptions);
    
	for vertex = 1:nVertices

		% Determine which vertexes are within the radius of the currently-picked vertex
        
		verticesCurrentlyWithinRadius = adjacencyMatrix(vertex,:);
		verticesCurrentlyWithinRadius = verticesCurrentlyWithinRadius(~isnan(verticesCurrentlyWithinRadius));        
        verticesCurrentlyWithinRadius = [vertex, verticesCurrentlyWithinRadius]; % add by Li Su 1-02-2012
        
		%n(vertex) = numel(verticesCurrentlyWithinRadius);

		for t = 1:nTimePoints

			% Work out the current time window
			currentTimeStart = (t - 1) * userOptions.temporalSearchlightResolution + 1;
			currentTimeWindow = (currentTimeStart : currentTimeStart + userOptions.temporalSearchlightWidth - 1);

			% Restrict data to current searchlight
			currentData = singleMesh(verticesCurrentlyWithinRadius, currentTimeWindow, :, :); % (vertices, time, condition, session)

			% Median over the time window
			switch lower(userOptions.searchlightPatterns)
				case 'spatial'
					% Spatial patterns: median over time window
					currentData = median(currentData, 2); % (vertices, 1, conditions, sessions)
					currentData = squeeze(currentData); % (vertices, conditions, sessions);
				case 'temporal'
					% Temporal patterns: mean over vertices within searchlight
					currentData = mean(currentData, 1); % (1, timePoints, conditions, sessions)
					currentData = squeeze(currentData); % (timePionts, conditions, sessions)
				case 'spatiotemporal'
					% Spatiotemporal patterns: all the data concatenated
					currentData = reshape(currentData, [], size(currentData, 3), size(currentData, 4)); % (dataPoints, conditions, sessions)
			end%switch:userOptions.sensorSearchlightPatterns

			searchlightRDM = zeros(nConditions, nConditions);

			% Average across sessions
			
			for session = 1:userOptions.nSessions

				searchlightRDM = searchlightRDM + squareform(pdist(squeeze(currentData(:,:,:,session))',userOptions.distance));
			
			end%for:sessions
			
			searchlightRDM = searchlightRDM / userOptions.nSessions;
		
			searchlightRDM = vectorizeRDM(searchlightRDM);
			
			% Locally store the full brain's worth of indexed RDMs.
			% searchlightRDMs(:,:,vertex, t) = squareform(searchlightRDM);
		
			try
				[rs, ps] = corr(searchlightRDM', modelRDMs_utv', 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
			catch
				[rs, ps] = corr(searchlightRDM', modelRDMs_utv, 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
			end%try
		


			smm_ps(vertex, t, :) = ps;
			smm_rs(vertex, t, :) = rs;

		end%for:t
		
		% Indicate progress every once in a while...
		if mod(vertex, floor(nVertices/20)) == 0, fprintf('.'); end%if
		
	end%for:sensorSite

    if userOptions.fisher
        smm_rs = fisherTransform(smm_rs);
    end%if
    
end%function