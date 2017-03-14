% Recipe_fMRI_searchlight
%
% Cai Wingfield 11-2009, 2-2010, 3-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = 'toolboxPathOnYourMachine'; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code
userOptions = defineUserOptions();

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%

fullBrainVols = rsa.fmri.fMRIDataPreparation('SPM', userOptions);
binaryMasks_nS = rsa.fmri.fMRIMaskPreparation(userOptions);

%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

models = rsa.constructModelRDMs(modelRDMs(), userOptions);

%%%%%%%%%%%%%%%%%
%% Searchlight %%
%%%%%%%%%%%%%%%%%

% prepare the searchlight RDMs
rsa.fmri.fMRIPrepareSearchlightRDMs(fullBrainVols, binaryMasks_nS, userOptions);

rsa.fmri.fMRISearchlightModelComparison(models, 'SPM', userOptions);
