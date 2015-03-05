function [smm_rs, smm_ps, searchlightRDMs] = searchlightMapping_MEG_source(singleMesh, Models, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%  [smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_MEG_source(singleMesh, Models, mask, userOptions, localOptions)
%  based on Li Su's script
% CW 5-2010, last updated by Li Su 3-2012

%% Get parameters

% Prepare Models
nModels = size(Models, 2);

if nModels>1
    modelRDMs_utv = squeeze(unwrapRDMs(vectorizeRDMs(Models)))'; % note by IZ: if Models has 1 model, this is a column vector, else models are put into rows
else
    modelRDMs_utv = squeeze(unwrapRDMs(vectorizeRDMs(Models)));
end


if userOptions.partial_correlation
    control_for_modelRDMs = modelRDMs_utv(userOptions.partial_modelNumber{1}, :);
    for m = 2:size(userOptions.partial_modelNumber, 2)
        control_for_modelRDMs = [control_for_modelRDMs;modelRDMs_utv(userOptions.partial_modelNumber{m}, :)];
    end
end

% Picking the right model
modelNumber = userOptions.modelNumber;
if nModels > 1 % added by IZ 11-12 to ensure model has more than 1 models to pick from
    modelRDMs_utv = modelRDMs_utv(modelNumber, :); % added by Li Su 11-2012 in order to test a single model at a time. not best in performance but keep thing simple.
end

nVertices = userOptions.nVertices;
nConditions = userOptions.nConditions;
%nSessions = userOptions.nSessions;

% How long is the stimulus (in data points)?
epochLength = size(singleMesh, 2); % (vertex, TIME, condition, session)

% Number of DATA points to loop with given the width and time step of
% searchlight updated by IZ 09-12
nTimePoints = floor((epochLength - (userOptions.temporalSearchlightWidth * userOptions.toDataPoints)) / ...
    (userOptions.temporalSearchlightResolution * userOptions.toDataPoints * userOptions.temporalDownsampleRate));


%% similarity-graph-map the volume with the searchlight

%n = nan(nVertices);
%searchlightRDMs = nan([nConditions, nConditions, nVertices, nTimePoints]);

% vertices change on the basis of maskng flag's value IZ 11-12
% updated: all searchlight run as masks IZ 03/12
vertices = userOptions.maskIndices.(userOptions.chi);


% Preallocate looped matrices for speed
smm_ps = zeros([nVertices, nTimePoints, nModels]);
smm_rs = zeros([nVertices, nTimePoints, nModels]);

% if strcmp(userOptions.groupStats,'FFX')
%     searchlightRDMs = single(zeros(nConditions,nConditions, nVertices, nTimePoints));
% end

for k = 1:length(vertices)
    vertex = vertices(k);
    % Determine which vertexes are within the radius of the currently-picked vertex
    
    verticesCurrentlyWithinRadius = userOptions.adjacencyMatrix(vertex,:);
    
    % remove nans
    verticesCurrentlyWithinRadius = verticesCurrentlyWithinRadius(~isnan(verticesCurrentlyWithinRadius));
    
    % add current vertex
    verticesCurrentlyWithinRadius = [vertex, verticesCurrentlyWithinRadius]; % add by Li Su 1-02-201
    
    % If masks are used, finding corresponding mask indices - update IZ 11-12
    if userOptions.maskingFlag
        location = ismember(verticesCurrentlyWithinRadius, vertices);
        verticesCurrentlyWithinRadius= verticesCurrentlyWithinRadius(location);
    end
    
    for t = 1:nTimePoints+1
        
        % Work out the current time window
        % converted to data points - updated by IZ 09-12
        currentTimeStart = (t - 1) * ...
            (userOptions.temporalSearchlightResolution * userOptions.temporalDownsampleRate * userOptions.toDataPoints) + 1;
        currentTimeWindow = ceil((currentTimeStart : currentTimeStart + ...
            (userOptions.temporalSearchlightWidth * userOptions.toDataPoints) - 1));
        
        currentData = singleMesh(verticesCurrentlyWithinRadius, currentTimeWindow, :, :); % (vertices, time, condition, session)
        
        searchlightRDM = zeros(nConditions, nConditions);
        
        % Average across sessions
        
        if not(userOptions.regularized)
            
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
            
            for session = 1:userOptions.nSessions
                
                searchlightRDM = searchlightRDM + squareform(pdist(squeeze(currentData(:,:,session))',userOptions.distance));
                
            end%for:sessions
            
        else % data regularization based on algorithm by Diedrichson et al 2011 - updated 12-12 IZ
            tempMesh = reshape(currentData, [], size(currentData, 3), size(currentData, 4));
            currentData = zeros(size(tempMesh, 1), size(tempMesh, 2) * size(tempMesh, 3)); % (data, conditions, sessions)
            
            % combining session-wise trials
            kk = 1;
            for j = 1:size(tempMesh,2)
                for i = 1:userOptions.nSessions
                    currentData(:, kk) = (tempMesh(:, j, i));
                    kk = kk + 1;
                end
            end
            
            r_matrix = g_matrix(zscore(squeeze(currentData(:,:)))', userOptions.nConditions, size(currentTimeWindow,2));
            searchlightRDM = searchlightRDM + (1 - r_matrix);
            
            if isnan(searchlightRDM) % sessions and conditions should be optimal
                error('Cannot calculate covariance matrix. Try reducing number of conditions');
            end
        end
        
        searchlightRDM = searchlightRDM / userOptions.nSessions;
        
        searchlightRDM = vectorizeRDM(searchlightRDM);
        
        % Locally store the full brain's worth of indexed RDMs. (just
        % lower triangle) added by IZ 09-12
        if strcmp(userOptions.groupStats, 'FFX')
%             searchlightRDMs(:,:,vertex, t) = single(tril(squareform(searchlightRDM)));
              searchlightRDMs.(['v_' num2str(vertex)]).(['t_' num2str(t)]).RDM = searchlightRDM;
        else
            searchlightRDMs = nan(1);
        end
        
        if userOptions.partial_correlation
            [rs, ps] = partialcorr(searchlightRDM', modelRDMs_utv', control_for_modelRDMs', 'type', userOptions.distanceMeasure, 'rows','pairwise');
        else
            try
                [rs, ps] = corr(searchlightRDM', modelRDMs_utv', 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
            catch ex
                [rs, ps] = corr(searchlightRDM', modelRDMs_utv, 'type', userOptions.distanceMeasure, 'rows', 'pairwise');
            end%try
        end
        
        smm_ps(vertex, t, :) = ps;
        smm_rs(vertex, t, :) = rs;
        
    end%for:t
    
    % Indicate progress every once in a while...
    if mod(vertex, floor(length(vertices) / 20)) == 0, fprintf('.'); end%if
    
end%for:vertices

if userOptions.fisher
    smm_rs = fisherTransform(smm_rs);
end%if

end%function
