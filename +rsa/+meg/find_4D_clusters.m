% Label spatio-temporal clusters in statistical maps.
% Created by Li Su, last update 17 Jun 2011

function [max_cluster_mass,lh_max_cluster_mass,rh_max_cluster_mass] = find_4D_clusters(inputFileName,outputFileName,connectivity_matrix,indexMasks,overwriteFlag,vertex_level_threshold)


import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% Threshold the statistical map, significant vertex lablled as 1 and non-sig vertex
%% lablled as 0. This results in activated vertex at the first level.

%modelName = [masked,'perm-r_value-',num2str(perm+range(1)-1)];
lh_inputFileName = [inputFileName,'-lh.stc']; % Generate filename
rh_inputFileName = [inputFileName,'-rh.stc']; % Generate filename  
lh_Vol = mne_read_stc_file1(lh_inputFileName); % Pull in data, requires MNE in the search path
rh_Vol = mne_read_stc_file1(rh_inputFileName); % Pull in data, requires MNE in the search path
lh_outputFileName = [outputFileName, '_cluster-lh.stc'];
rh_outputFileName = [outputFileName, '_cluster-rh.stc'];

% do we need this check here? Commenting it out --IZ 11-12
% if (exist(lh_outputFileName, 'file') && exist(rh_outputFileName, 'file')) %&& not(overwriteFlag)
%     fprintf(['\b\b cluster files (', lh_outputFileName, ') already exist', '. Skip... ']);
%     lh_cluster = mne_read_stc_file(lh_outputFileName);
%     rh_cluster = mne_read_stc_file(rh_outputFileName);
    
% else
    time_window = 1:size(lh_Vol.data,2);
    % remove 1s in the end of the matrix and change back to p from 1-p. may not need this two lines in the final version.
    lh = lh_Vol.data(:,time_window);
    rh = rh_Vol.data(:,time_window);
    
    % commenting this out as the masking is alrady done in previous stages
    % IZ 03/12
%     if size(indexMasks) > 0
%         disp('applying masking...');
%         nMasks = numel(fieldnames(indexMasks));
%         lh = MEG_masking(lh,indexMasks,1:nMasks); % I dont think this works fine without prior chirality info 03/12 IZ
%         rh = MEG_masking(rh,indexMasks,1:nMasks);
%     end

    lh(lh < vertex_level_threshold) = 0;
    lh(lh >= vertex_level_threshold) = 1;
    rh(rh < vertex_level_threshold) = 0;
    rh(rh >= vertex_level_threshold) = 1;


    %% find cluster at each time point    
    lh_cluster = lh_Vol;
    rh_cluster = rh_Vol;
    lh_cluster.data = [];
    rh_cluster.data = [];

    for t = time_window
        if t==floor(size(time_window,2)/2)
            fprintf('.');
        end
        lh_masked_connectivity_matrix = connectivity_matrix;
        rh_masked_connectivity_matrix = connectivity_matrix;

        % select only the activated vertex from connectivity matrix and
        % disconnect all inactive vertex.
        
        % Following implementation uses the connectivity matrix as a sparse 
        % matrix to perform assignments more efficiently
        lh_masked_connectivity_matrix(lh(:,t) == 0,:) = 0;
        lh_masked_connectivity_matrix(:,lh(:,t) == 0) = 0;
        rh_masked_connectivity_matrix(rh(:,t) == 0,:) = 0;
        rh_masked_connectivity_matrix(:,rh(:,t) == 0) = 0;
        
        % converting back to full matrix
        lh_masked_connectivity_matrix = full(lh_masked_connectivity_matrix);
        rh_masked_connectivity_matrix = full(rh_masked_connectivity_matrix);
        
        % alternatively, previous step can be performed on full matrix as well
%         lh_masked_connectivity_matrix(lh(:,t) == 0,:) = 0;
%         lh_masked_connectivity_matrix(:,lh(:,t) == 0) = 0;
%         rh_masked_connectivity_matrix(rh(:,t) == 0,:) = 0;
%         rh_masked_connectivity_matrix(:,rh(:,t) == 0) = 0;
                
        [lh_p,lh_q,lh_r,lh_s] = dmperm(lh_masked_connectivity_matrix);
        [rh_p,rh_q,rh_r,rh_s] = dmperm(rh_masked_connectivity_matrix);

        clear lh_s rh_s;

        z_lh = zeros(size(connectivity_matrix,1),1);
        z_rh = zeros(size(connectivity_matrix,1),1);

        for i = 1:size(lh_r,2)-2
            cluster = lh_p(lh_r(i):lh_r(i+1)-1);            
            if t > 1 % check if this cluster overlaps with the previous time point
                linked_cluster = intersect(find(lh_cluster.data(:,t-1) ~= 0),cluster); 
                if ~isempty(linked_cluster)
                    backwards_cluster_IDs = lh_cluster.data(linked_cluster,t-1);
                    current_cluster_ID = min(backwards_cluster_IDs);
                    lh_cluster.data(ismember(lh_cluster.data(:,1:t-1),backwards_cluster_IDs)) = current_cluster_ID;
                    z_lh(cluster) = current_cluster_ID;
                else
                    z_lh(cluster) = max(max(max(max(lh_cluster.data(:,1:t-1))),z_lh))+1;
                end
            else
                z_lh(cluster) = i;
            end
        end

        for i = 1:size(rh_r,2)-2
            cluster = rh_p(rh_r(i):rh_r(i+1)-1);            
            if t > 1 % check if this cluster overlaps with the previous time point
                linked_cluster = intersect(find(rh_cluster.data(:,t-1) ~= 0),cluster); 
                if ~isempty(linked_cluster)
                    backwards_cluster_IDs = rh_cluster.data(linked_cluster,t-1);
                    current_cluster_ID = min(backwards_cluster_IDs);
                    rh_cluster.data(ismember(rh_cluster.data(:,1:t-1),backwards_cluster_IDs)) = current_cluster_ID;
                    z_rh(cluster) = current_cluster_ID;
                else
                    z_rh(cluster) = max(max(max(max(rh_cluster.data(:,1:t-1))),z_rh))+1;
                end
            else
                z_rh(cluster) = i;
            end
        end
        lh_cluster.data(:,t) = z_lh;
        rh_cluster.data(:,t) = z_rh;
        
        % Indicate progress every once in a while...
        if mod(t, floor(max(time_window)/50)) == 0, fprintf('\b.:'); end%if

    end
    
    mne_write_stc_file1(lh_outputFileName, lh_cluster);
    mne_write_stc_file1(rh_outputFileName, rh_cluster);

%% compute the cluster mass
lh_data = reshape(lh_cluster.data,1,size(lh_cluster.data,1)*size(lh_cluster.data,2));
lh_t_values = reshape(lh_Vol.data,1,size(lh_Vol.data,1)*size(lh_Vol.data,2));
numberOfCluster_left = max(lh_data); 
lh_max_cluster_mass = 0;
for i = 1:numberOfCluster_left 
    cluster_mass = sum(lh_t_values(lh_data == i)-vertex_level_threshold);
    if cluster_mass > lh_max_cluster_mass
        lh_max_cluster_mass = cluster_mass;
    end
end

rh_data = reshape(rh_cluster.data,1,size(rh_cluster.data,1)*size(rh_cluster.data,2));
rh_t_values = reshape(rh_Vol.data,1,size(rh_Vol.data,1)*size(rh_Vol.data,2));
numberOfCluster_right = max(rh_data);
rh_max_cluster_mass = 0;
for i = 1:numberOfCluster_right 
    cluster_mass = sum(rh_t_values(rh_data == i)-vertex_level_threshold);
    if cluster_mass > rh_max_cluster_mass
        rh_max_cluster_mass = cluster_mass;
    end
end

max_cluster_mass = max(lh_max_cluster_mass,rh_max_cluster_mass);    
