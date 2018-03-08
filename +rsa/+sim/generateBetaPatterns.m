function OUT = generateBetaPatterns(clusterSpec, nVoxels, centroid)

% B = generateBetaPatterns(clusterSpec, nVoxels[, centroid])
%
% Generates a [nConditions nVoxels]-sized beta matrix of data where the
% responses are clustered according to clusterSpec (see below). Does not
% simulate scanner noise --- these simulated responses are "true".
%
%        B --- The beta matrix.
%                A [nConditions nVoxels]-sized matrix of simulated response.
%
%        clusterSpec --- Specification for the clustering of conditions.
%                clusterSpec is a cell array. Clustering of the conditions is
%                done in the following manner:
%                - Each cell represents a hierarchy.
%                - The first entry in a cell is the 'spread' of that level of
%                  the hierarchy.
%                - If the second entry is another cell, it and all following
%                  entries (which must also be cells) are sub-hierarchies.
%                - If the second entry is a number (in which case there must
%                  only  be two entries), this number represents a number of
%                  'leaves', such that the cell is hierarchically an 'atom'.
%                So we have a well-defined recursive datatype for hierarchies.
%                For example,
%                        clusterSpec = {20, {6, {3, 5},{2, 3}}, {4, 7}}
%                represents the following hierarchy:
%                                            |
%                                ------------------------- 20
%                                |                       |
%                          ------------- 6               |
%                          |           |              ------- 4
%                        ----- 3      --- 2           |||||||
%                        |||||        |||             |||||||
%                         [5]         [3]               [7]
%
%        nVoxels --- The number of voxels in the simulated RoI.
%                An integer.
%
%        centroid --- The centre of the clusters.
%                This is mainly used for the recursion in this function, but
%                a vector of length nVoxels put here will centre the clusters
%                around this point in voxel space. Defaults to the origin.
%
% Cai Wingfield 6-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

	% Probably on the first call, there is no centroid specified (an
	% implicit origin). Later, centroids will be picked based on the
	% results of previous iterations.
	if ~exist('centroid', 'var')
		centroid = zeros(1, nVoxels);
	end%if
	
	% The first is the spread for this level of the hierarchy.
	sigma = sqrt(clusterSpec{1});
	
	% If there are no subclusters...
	if ~iscell(clusterSpec{2})
		
		% The number of conditions at this level is the number of
		% NaN-filled rows.
		nConditions = clusterSpec{2};
		
		% These conditions are given points at a spread of sigma around the
		% current centroid.
		allOnCentroid = repmat(centroid, nConditions, 1);
		OUT = allOnCentroid + sigma * randn(nConditions, nVoxels);
		
	% Otherwise, there are subclusters
	else
		
		% There are as many subclusters as elements of this cell (minus
		% 1, for the sigma).
		nSubclusters = numel(clusterSpec) - 1; 
		
		OUT = [];
		
		% For each subcluster...
		for subclusterIndex = 1:nSubclusters
			
			% The centroid for this subcluster will be the centroid of the
			% parent cluster plus a random pertubation.
			subclusterCentroid = centroid + sigma * randn(1, nVoxels);
			
			% All further subclusters are now recursively generated,
			% centred about this local centroid.
			OUT = [OUT; generateBetaPatterns(clusterSpec{subclusterIndex + 1}, nVoxels, subclusterCentroid)];
			
		end%for:subcluster
		
	end%if:top-left only
	
end%function