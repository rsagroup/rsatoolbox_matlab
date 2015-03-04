function [max_clust_mass, clust_mass]=find_max_clust_mass(data,thresh)

%% Function called by permutation_cluster_test_2dtfr_func.m
%
% data is matrix of t-values, thresh is t-value (or other test statistic)
%
% Written by A.Ghuman (2009), adapted by Alex C (07/09)

%% Finds elements in data above threshold
bwmatrix = zeros(size(data));
bwmatrix(find(data>thresh)) = 1;

% Numbers clusters
clusters_matrix = bwlabel(bwmatrix);

% Determines extent of each cluster
for i = 1:max(max(clusters_matrix))
   clust_mass(i) = sum(data(find(clusters_matrix==i)));
end

% Only takes max cluster
if ~isempty(find(bwmatrix==1))
    max_clust_mass=max(clust_mass);
else
    clust_mass=0;
    max_clust_mass=0;
end
