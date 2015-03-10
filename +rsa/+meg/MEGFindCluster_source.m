% Label spatio-temporal clusters in statistical maps.
% Created by Li Su, last update 15-02-2012


function MEGFindCluster_source(Models, range, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

% The analysisName will be used to label the files which are eventually saved.
modelNumber = userOptions.modelNumber;
modelName = spacesToUnderscores(Models(modelNumber).name);
if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end
MapsFilename = ['perm-', num2str(max(range)), '_', modelName, '_t_map_cluster'];
if userOptions.maskingFlag
    MapsFilename = [MapsFilename , '_masked'];
end

if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end

promptOptions.functionCaller = 'MEGFindCluster_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-lh.stc']);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-rh.stc']);

% TODO: This prevents parfor toolbox from working
overwriteFlag = true;%overwritePrompt(userOptions, promptOptions);

if overwriteFlag
        
    inputFileName = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_tMesh_' modelName '_allSubjects']);
    if userOptions.maskingFlag
        inputFileName = [inputFileName , '_masked'];
    end
    
    if isnan(userOptions.primaryThreshold)
        % find the vertex level threshold of t value by sorting the t
        % distribution and select the cut off point of x%
        % percent = userOptions.primaryThreshold;
        
        t_values_lh = mne_read_stc_file1([inputFileName,'-lh.stc']);
        t_values_rh = mne_read_stc_file1([inputFileName,'-rh.stc']);
        
        t_values_lh.data = reshape(t_values_lh.data,1,size(t_values_lh.data,1)*size(t_values_lh.data,2));
        t_values_rh.data = reshape(t_values_rh.data,1,size(t_values_rh.data,1)*size(t_values_rh.data,2));
        
        if userOptions.maskingFlag
            nMasks = size(userOptions.maskNames,2);
            t_values_lh.data = MEG_masking(t_values_lh.data, userOptions.indexMasks, 1:2:nMasks); % assuming oddly numbered are left hemisphere and vice versa IZ 03/12
            t_values_rh.data = MEG_masking(t_values_rh.data, userOptions.indexMasks, 2:2:nMasks);
        end
        
        t_distribution = sort([t_values_lh.data,t_values_rh.data]);
        percent = 0.05;
        vertex_level_threshold = t_distribution(ceil(size(t_distribution,2)*(1-percent)));
    else
        %% set the threshold to p < .05 at the vertex level according to
        %% random effect T-distribution (one tailed). df = number of
        %% subjects - 1.
        pval = userOptions.primaryThreshold;
        vertex_level_threshold = abs(tinv(pval,(userOptions.nSubjects-1)));
    end
    
    %% Build a connectivity (sparse) matrix representing the adjacency
    %% informaiton between vertex.
    
    number_of_vertex = userOptions.targetResolution;
    minDist = userOptions.minDist;
    
    adjMatrix = calculateMeshAdjacency(number_of_vertex, minDist, userOptions);
    
    start_vertex = repmat(1:number_of_vertex,1,7)';
    end_vertex = vertcat(reshape(adjMatrix, number_of_vertex*6,1), (1:number_of_vertex)');
    adjVect = [start_vertex, end_vertex];
    adjVect = adjVect(~any(isnan(adjVect),2),:);
    connectivity_matrix = sparse(adjVect(:,1), adjVect(:,2), 1);
    % convert to full or sparse depending on the script of
    % find_4D_clusters(...)
    % connectivity_matrix = full(connectivity_matrix);
    % connectivity_matrix = connectivity_matrix - eye(number_of_vertex);
    
    fprintf('Performing 4D spatiotemporal clustering...\n');
    fprintf('.:');
    
    find_4D_clusters(inputFileName,inputFileName,connectivity_matrix,userOptions.indexMasks, overwriteFlag, vertex_level_threshold);
    
    %range = 1:userOptions.significanceTestPermutations; % in the future, the follows should be run in multiple parallel CPUs
    
    numberOfPermutation = size(range,2);
    
    for perm = min(range):max(range)
        % disp([num2str(round(perm/numberOfPermutation*100)) '% finished.  ']);
        fprintf('.');
        inputFileName = fullfile(userOptions.rootPath, 'Maps', modelName, ['perm-' num2str(perm) '_' modelName '_t_map']);
        
        [max_cluster_mass(perm-min(range)+1),lh_max_cluster_mass(perm-min(range)+1),rh_max_cluster_mass(perm-min(range)+1)] ...
            = find_4D_clusters(inputFileName,inputFileName,connectivity_matrix,userOptions.indexMasks, overwriteFlag, vertex_level_threshold);
        
        if ~userOptions.debug
            delete([inputFileName,'-lh.stc']);
            delete([inputFileName,'-rh.stc']);
        end
        
    end
    
    csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution' num2str(min(range)) '.csv']),max_cluster_mass);
    %csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution-lh' num2str(min(range)) '.csv']),lh_max_cluster_mass);
    %csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution-rh' num2str(min(range)) '.csv']),rh_max_cluster_mass);
    fprintf(' null distribution saved.\n');
    
    %    timeStamp = datestr(now);
    %     gotoDir(userOptions.rootPath, 'Details');
    %     cd(fullfile(userOptions.rootPath, 'Details'));
    %     fprintf(['Saving Details to ' fullfile(pwd, DetailsFilename) '\n']);
    %     save(DetailsFilename, 'timeStamp', 'userOptions');
    
    %     else % masking flag set to true IZ 12-12
    %
    %         inputFileName = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_tMesh_' modelName '_allSubjects_masked']);
    %
    %         if isnan(userOptions.primaryThreshold)
    %             % find the vertex level threshold of t value by sorting the t
    %             % distribution and select the cut off point of x%
    %             % percent = userOptions.primaryThreshold;
    %
    %             t_values_lh = mne_read_stc_file1([inputFileName,'-lh.stc']);
    %             t_values_rh = mne_read_stc_file1([inputFileName,'-rh.stc']);
    %
    %             t_values_lh.data = reshape(t_values_lh.data,1,size(t_values_lh.data,1)*size(t_values_lh.data,2));
    %             t_values_rh.data = reshape(t_values_rh.data,1,size(t_values_rh.data,1)*size(t_values_rh.data,2));
    %
    %             %                 if size(indexMasks) > 0
    %             %                     nMasks = size(userOptions.maskNames,2);
    %             %                     t_values_lh.data = MEG_masking(t_values_lh.data, indexMasks, 1:nMasks);
    %             %                     t_values_rh.data = MEG_masking(t_values_rh.data, indexMasks, 1:nMasks);
    %             %                 end
    %
    %             t_distribution = sort([t_values_lh.data,t_values_rh.data]);
    %             percent = 0.05;
    %             vertex_level_threshold = t_distribution(ceil(size(t_distribution,2)*(1-percent)));
    %         else
    %             %% set the threshold to p < .05 at the vertex level according to
    %             %% random effect T-distribution (one tailed). df = number of
    %             %% subjects - 1.
    %             pval = userOptions.primaryThreshold;
    %             vertex_level_threshold = abs(tinv(pval,(userOptions.nSubjects-1)));
    %         end
    %
    %         %% Build a connectivity (sparse) matrix representing the adjacency
    %         %% informaiton between vertex.
    %
    %         number_of_vertex = userOptions.targetResolution;
    %         minDist = userOptions.minDist;
    %
    %         adjMatrix = calculateMeshAdjacency(number_of_vertex, minDist, userOptions);
    %
    %         start_vertex = repmat(1:number_of_vertex,1,7)';
    %         end_vertex = vertcat(reshape(adjMatrix, number_of_vertex*6,1), (1:number_of_vertex)');
    %         adjVect = [start_vertex, end_vertex];
    %         adjVect = adjVect(~any(isnan(adjVect),2),:);
    %
    %         %             k=1;
    %         %             for i=1:size(indexMasks.bankssts_lh.tw_101_200.maskIndices,1)
    %         %                 x = find(adjVect(:,1)==indexMasks.(thisMask_lh).(thisTimeWindow).maskIndices(i));
    %         %                 for j =1:size(x,1)
    %         %                     if adjVect(x(j),2)==indexMasks.(thisMask_lh).(thisTimeWindow).maskIndices
    %         %                         adjVect_mask(k,:) = adjVect(x,:);
    %         %                         k=k+1;
    %         %                     end
    %         %                 end
    %         %                 adjVect(
    %         %             end
    %         connectivity_matrix = sparse(adjVect(:,1), adjVect(:,2), 1);
    %         % convert to full or sparse depending on the script of
    %         % find_4D_clusters(...)
    %         % connectivity_matrix = full(connectivity_matrix);
    %         % connectivity_matrix = connectivity_matrix - eye(number_of_vertex);
    %
    %         fprintf('Performing 4D spatiotemporal clustering...\n');
    %         fprintf('.:');
    %
    %         find_4D_clusters(inputFileName,inputFileName,connectivity_matrix,currentMasks, overwriteFlag, vertex_level_threshold);
    %
    %         for perm = min(range):max(range)
    %             fprintf('.');
    %             inputFileName = fullfile(userOptions.rootPath, 'Maps', modelName, ['perm-' num2str(perm) '_' modelName '_t_map']);
    %
    %             [max_cluster_mass(perm-min(range)+1),lh_max_cluster_mass(perm-min(range)+1),rh_max_cluster_mass(perm-min(range)+1)] ...
    %                 = find_4D_clusters(inputFileName,inputFileName,connectivity_matrix,currentMasks, overwriteFlag, vertex_level_threshold);
    %
    %             if ~userOptions.debug
    %                 delete([inputFileName,'-lh.stc']);
    %                 delete([inputFileName,'-rh.stc']);
    %             end
    %
    %         end
    %
    %         csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution' num2str(min(range)) '.csv']),max_cluster_mass);
    %         %csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution-lh' num2str(min(range)) '.csv']),lh_max_cluster_mass);
    %         %csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution-rh' num2str(min(range)) '.csv']),rh_max_cluster_mass);
    %         fprintf(' null distribution saved.\n');
    %
    %         %    timeStamp = datestr(now);
    %         %     gotoDir(userOptions.rootPath, 'Details');
    %         %     cd(fullfile(userOptions.rootPath, 'Details'));
    %         %     fprintf(['Saving Details to ' fullfile(pwd, DetailsFilename) '\n']);
    %         %     save(DetailsFilename, 'timeStamp', 'userOptions');
    %         %end % for mask
    %     end % if maskingFlag
else
    
end
cd(returnHere); % And go back to where you started
