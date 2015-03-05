% Recipe_MEG_searchlight_source (which_model)
%
% which_model: the model number, in the order of they occur in the modelRDMs.m
%
% Cai Wingfield 5-2010, 8-2010
% update by Li Su 3-2012, 11-2012
% updated Fawad 12-2013, 02-2014, 10-2014

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
toolboxRoot = 'C:\Users\cai\code\rsagroup-rsatoolbox\'; 
addpath(genpath(toolboxRoot));
userOptions = defineUserOptions();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
userOptions.modelNumber = which_model;
models = rsa.constructModelRDMs(userOptions);

%%%%%%%%%%%%%%%%%%
%% Set metadata %%
%%%%%%%%%%%%%%%%%%
% TODO: This is bad
userOptions = rsa.meg.setMetadata_MEG(models, userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Mask preparation %% 
%%%%%%%%%%%%%%%%%%%%%%
% TODO: Why is this set in the userOptions struct?
if userOptions.maskingFlag
    userOptions.indexMasks = rsa.meg.MEGMaskPreparation_source(userOptions);  
else
    userOptions.indexMasks = rsa.meg.allBrainMask(userOptions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate adjacency matrix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Why is this stored in the userOptions struct separately to the list
% TODO: of subjects?
nSubject = userOptions.nSubjects;
userOptions.adjacencyMatrix = rsa.meg.calculateMeshAdjacency(userOptions.nVertices, userOptions.sourceSearchlightRadius, userOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Starting parallel toolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if userOptions.flush_Queue
    rsa.par.flushQ();
end

if userOptions.run_in_parallel
    rsa.par.initialise_CBU_Queue(userOptions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Searchlight - Brain RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Consider a timestamped printing function.
parfor subject = 1:nSubject
    rsa.meg.MEGSearchlight_source(subject, models, userOptions);
end
fprintf('Stage 1 - Searchlight Brain RDM Calculation: ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smoothing and upsampling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: What is `stages`?
%if stages == 2
%    interpolate_MEG_source(Models, userOptions);
%end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Random permutation %%
%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(userOptions.groupStats, 'FFX')
    tic
    % fixed effect test
    fprintf('Averaging RDMs across subjects and performing permutation tests to calculate r-values.');
    rsa.meg.FFX_permutation(models, userOptions)
    fprintf('Stage 2 - Fixed Effects Analysis: ');
    toc
else
    % random effect test
    fprintf('Performing permutation tests over all subjects (random effect) to calculate p-values.');
    tic
    rsa.meg.RFX_permutation(models,userOptions)
    fprintf('Stage 2 - Random Effects Analysis: ');
    toc 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sptiotemporal clustering %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jobSize = userOptions.jobSize;
number_of_permutations = userOptions.significanceTestPermutations;

tic
parfor j = 1:number_of_permutations/jobSize
    range = (j-1)*jobSize+1:j*jobSize;
    rsa.meg.MEGFindCluster_source(models, range, userOptions);
end
fprintf('Stage 3 - Spatiotemporal Clustering: ' );
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute cluster level p-values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
rsa.meg.get_cluster_p_value(models,userOptions);
fprintf('Stage 4 - Computing cluster level p values: ');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stopping parallel toolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if userOptions.run_in_parallel
    matlabpool close;
end

%%%%%%%%%%%%%%%%%%%%
%% Delete Selected Directories%%
%%%%%%%%%%%%%%%%%%%%
if (userOptions.deleteTMaps_Dir || userOptions.deleteImageData_Dir || userOptions.deletePerm)
    rsa.util.deleteDir(userOptions, models);
end

%%%%%%%%%%%%%%%%%%%%
%% Sending an email %%
%%%%%%%%%%%%%%%%%%%%
if userOptions.recieveEmail
    rsa.par.setupInternet();
    rsa.par.setupEmail(userOptions.mailto);
end
