%  Recipe_MEG_sensor
%
% Cai Wingfield 12-2009, 8-2010
% Updated by IZ 07/13

function Recipe_MEG_sensor(which_model)

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = '/imaging/iz01/test/toolbox/devel/toolbox/'; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code
userOptions = projectOptions();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
userOptions.modelNumber = which_model;
Models = constructModelRDMs(userOptions);

%%%%%%%%%%%%%%%%%%%%
%% Set Meta-data %%
%%%%%%%%%%%%%%%%%%%%
userOptions = setMetadata_MEG(Models, userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%

sensorImages = MEGDataPreparation_sensor(userOptions);
maskedSensors = MEGDataMasking_sensor(sensorImages, userOptions);


%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

RDMs = constructRDMs(maskedSensors, betaCorrespondence(), userOptions);
RDMs = averageRDMs_subjectSession(RDMs, 'session');
aRDMs = averageRDMs_subjectSession(RDMs, 'subject');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figureRDMs(aRDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1)); % Display the calculated RDMs
figureRDMs(Models, userOptions, struct('fileName', 'ModelRDMs', 'figureNumber', 2)); % Display the models

MDSConditions(aRDMs, userOptions);
dendrogramConditions(aRDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second-order analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

pairwiseCorrelateRDMs({aRDMs, Models}, userOptions);
MDSRDMs({aRDMs, Models}, userOptions);
distanceBarRDMs({RDMs}, {Models}, userOptions);

testSignificance({aRDMs}, {Models}, userOptions); % fixed effects analysis
testSignificance_RandomEffects(RDMs, Models, userOptions); % random effects analysis

end
