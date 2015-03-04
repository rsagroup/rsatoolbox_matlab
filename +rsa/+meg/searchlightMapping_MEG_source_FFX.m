function  searchlightRDMs = searchlightMapping_MEG_source_FFX (singleMesh, Models, userOptions)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

	%  [smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_MEG_source(singleMesh, Models, mask, userOptions, localOptions)
	%  
	%  Li Su 3-2012
    
	RDM_vector_size = max(size(unwrapRDMs(vectorizeRDMs(Models))));
    
	nVertices = userOptions.nVertices;
	nConditions = userOptions.nConditions;
	%nSessions = userOptions.nSessions;

	% How long is the stimulus (in time points)?
	epochLength = size(singleMesh, 2); % (vertex, time, condition, session)

	% Number of time points to loop with given the withd and time step of searchlight
	nTimePoints = floor((epochLength - userOptions.temporalSearchlightWidth) / userOptions.temporalSearchlightResolution);

	searchlightRDMs = nan([RDM_vector_size, nVertices, nTimePoints]);
    
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
			searchlightRDMs(:,vertex, t) = searchlightRDM;

		end%for:t
		
		% Indicate progress every once in a while...
		if mod(vertex, floor(nVertices/20)) == 0, fprintf('.'); end%if
		
	end%for:sensorSite
    
end%function
