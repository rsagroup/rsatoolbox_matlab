function simulationOptions = simulationOptions_demo_SL()
%
% simulationOptions_demo
%
% This function is used to make a simulationOptions struct for use in the
% simulation part of the demo.  The options should be set to personal
% preference.
%
% Cai Wingfield 7-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council
% A triple containing the dimensions of the simulated RoI in voxels.
simulationOptions.volumeSize_vox = [7 7 7];

% A triple containing the dimensions of one voxel in mm.
simulationOptions.voxelSize_mm = [3 3 3.75];

% The number of repititions of each stimulus which is "presented" to each of the
% simulated subjects.
simulationOptions.nRepititions = 3;% 10

% The duration of the stimulus presentation in seconds.
simulationOptions.stimulusDuration = 0.3;

% The duration of one trial in seconds.
simulationOptions.trialDuration = 3;

% The time for one TR in seconds.
simulationOptions.TR = 1.5; % seconds

% The amount of noise to be added by the simulated scanner. This corresponds to
% the square of the standard deviation of the gaussian distibution from which
% the noise is drawn (?).
simulationOptions.scannerNoiseLevel = 10000;% used to be 3000

% A 4-tuple. The first three entries are the x, y and z values for the gaussian
% spatial smoothing kernel FWHM in mm and the fourth is the size of the temporal
% smoothing FWHM.
simulationOptions.spatiotemporalSmoothingFWHM_mm_s = [4 4 4 4.5];

% The specification for the clustering of conditions.
% This should be a cell array. Clustering of the conditions is done in the
% following manner:
%        - Each cell represents a hierarchy.
%        - The first entry in a cell is the 'spread' of that level of the
%          hierarchy.
%        - If the second entry is another cell, it and all following entries
%          (which must also be cells) are sub-hierarchies.
%        - If the second entry is a number (in which case there must only be
%          two entries), this number represents a number of 'leaves', such that
%          the cell is hierarchically an 'atom'.
% So we have a well-defined recursive datatype for hierarchies.
% For example,
%        clusterSpec = {20, {6, {3, 5},{2, 3}}, {4, 7}}
% represents the following hierarchy:
%                                  |
%                      ------------------------- 20
%                      |                      |
%               ------------- 6               |
%               |           |              ------- 4
%             ----- 3      --- 2           |||||||
%             |||||        |||             |||||||
%               [5]         [3]               [7]
simulationOptions.clusterSpec = ...
	{2, ...
		{1, ...
			{.5, 10}, ...
			{.5, 10} ...
		}, ...
		{1, ...
			{.5, 10}, ...
			{.5, 10} ...
		} ...
	};
% for searchlight demo
simulationOptions.brainVol = [64 64 32];
simulationOptions.effectCen = [20 20 15];

simulationOptions.nConditions = 40;

end%function

