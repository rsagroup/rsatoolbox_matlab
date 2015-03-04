% Recipe_MEG_searchlight_sensor
%
% Cai Wingfield 3-2010, 8-2010

function Recipe_MEG_searchlight_sensor(which_model)

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
toolboxRoot = '/imaging/iz01/test/toolbox/devel/toolbox/'; addpath(genpath(toolboxRoot));
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

%%%%%%%%%%%%%%%%%
%% Searchlight %%
%%%%%%%%%%%%%%%%%
MEGSearchlight_sensor(sensorImages, Models, userOptions);

end
