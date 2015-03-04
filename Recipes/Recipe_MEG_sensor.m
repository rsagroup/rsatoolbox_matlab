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
Models = rsa.constructModelRDMs(userOptions);

%%%%%%%%%%%%%%%%%%%%
%% Set Meta-data %%
%%%%%%%%%%%%%%%%%%%%
userOptions = rsa.meg.setMetadata_MEG(Models, userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%

sensorImages  = rsa.meg.MEGDataPreparation_sensor(userOptions);
maskedSensors = rsa.meg.MEGDataMasking_sensor(sensorImages, userOptions);


%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

RDMs  = rsa.constructRDMs(maskedSensors, betaCorrespondence(), userOptions);
RDMs  = rsa.rdm.averageRDMs_subjectSession(RDMs, 'session');
aRDMs = rsa.rdm.averageRDMs_subjectSession(RDMs, 'subject');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order visualisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rsa.figureRDMs(aRDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1)); % Display the calculated RDMs
rsa.figureRDMs(Models, userOptions, struct('fileName', 'ModelRDMs', 'figureNumber', 2)); % Display the models

rsa.MDSConditions(aRDMs, userOptions);
rsa.dendrogramConditions(aRDMs, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second-order analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

rsa.pairwiseCorrelateRDMs({aRDMs, Models}, userOptions);
rsa.MDSRDMs({aRDMs, Models}, userOptions);
rsa.compareRefRDM2candRDMs(RDMs(1), Models, userOptions);

rsa.stat.testSignificance({aRDMs}, {Models}, userOptions); % fixed effects analysis
rsa.stat.testSignificance_RandomEffects(RDMs, Models, userOptions); % random effects analysis

end
