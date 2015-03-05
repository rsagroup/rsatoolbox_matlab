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
Models = rsa.constructModelRDMs(userOptions);

%%%%%%%%%%%%%%%%%%%%
%% Set Meta-data %%
%%%%%%%%%%%%%%%%%%%%
userOptions = rsa.meg.setMetadata_MEG(Models, userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%
sensorImages = rsa.meg.MEGDataPreparation_sensor(userOptions);

%%%%%%%%%%%%%%%%%
%% Searchlight %%
%%%%%%%%%%%%%%%%%
rsa.meg.MEGSearchlight_sensor(sensorImages, Models, userOptions);

end
