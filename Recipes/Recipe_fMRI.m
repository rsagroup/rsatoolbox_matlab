% Recipe_fMRI
% this 'recipe' performs region of interest analysis on fMRI data.
% Cai Wingfield 5-2010, 6-2010, 7-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = 'toolboxPathOnYourMachine'; addpath(genpath(toolboxRoot));
userOptions = defineUserOptions();

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%

fullBrainVols = rsa.core.fMRIDataPreparation('SPM', userOptions);
binaryMasks_nS = rsa.core.fMRIMaskPreparation(userOptions);
responsePatterns = rsa.core.fMRIDataMasking(fullBrainVols, binaryMasks_nS, 'SPM', userOptions);

%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

RDMs = rsa.constructRDMs(responsePatterns, 'SPM', userOptions);
sRDMs = rsa.core.averageRDMs_subjectSession(RDMs, 'session');
RDMs = rsa.core.averageRDMs_subjectSession(RDMs, 'session', 'subject');

Models = rsa.core.constructModelRDMs(modelRDMs(), userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rsa.figureRDMs(RDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1));
rsa.figureRDMs(Models, userOptions, struct('fileName', 'ModelRDMs', 'figureNumber', 2));

rsa.MDSConditions(RDMs, userOptions);
rsa.dendrogramConditions(RDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% relationship amongst multiple RDMs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rsa.pairwiseCorrelateRDMs({RDMs, Models}, userOptions);
rsa.MDSRDMs({RDMs, Models}, userOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistical inference %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
roiIndex = 1;% index of the ROI for which the group average RDM will serve 
% as the reference RDM. 
for i=1:numel(Models)
    models{i}=Models(i);
end
userOptions.RDMcorrelationType='Kendall_taua';
userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
userOptions.RDMrelatednessThreshold = 0.05;
userOptions.figureIndex = [10 11];
userOptions.RDMrelatednessMultipleTesting = 'FDR';
userOptions.candRDMdifferencesTest = 'subjectRFXsignedRank';
userOptions.candRDMdifferencesThreshold = 0.05;
userOptions.candRDMdifferencesMultipleTesting = 'none';
stats_p_r=rsa.compareRefRDM2candRDMs(RDMs(roiIndex), models, userOptions);
