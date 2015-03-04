% Recipe_MEG_searchlight_source (which_model)
%
% which_model: the model number, in the order of they occur in the modelRDMs.m
%
% Cai Wingfield 5-2010, 8-2010
% update by Li Su 3-2012, 11-2012
% updated Fawad 12-2013, 02-2014, 10-2014

function Recipe_MEG_searchlight_source (which_model)

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
toolboxRoot = '/home/fj01/toolbox'; 
addpath(genpath(toolboxRoot));
userOptions = projectOptions();

%%%%%%%%%%%%%%%%%%%%
%% Starting parallel toolbox %%
%%%%%%%%%%%%%%%%%%%%
if userOptions.flush_Queue
    flushQ();
end

if userOptions.run_in_parallel
    p = initialise_CBU_Queue(userOptions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
userOptions.modelNumber = which_model;
Models = constructModelRDMs(userOptions);

%%%%%%%%%%%%%%%%%%
%% Set metadata %%
%%%%%%%%%%%%%%%%%%
userOptions = setMetadata_MEG(Models, userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Mask preparation %% 
%%%%%%%%%%%%%%%%%%%%%%
if userOptions.maskingFlag
    userOptions.indexMasks = MEGMaskPreparation_source(userOptions);  
else
    userOptions.indexMasks = allBrainMask(userOptions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Searchlight - Brain RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSubject = userOptions.nSubjects;
tic
userOptions.adjacencyMatrix = calculateMeshAdjacency(userOptions.nVertices, userOptions.sourceSearchlightRadius, userOptions);
parfor subject = 1:nSubject
    MEGSearchlight_source(subject, Models, userOptions);
end
fprintf('Stage 1 - Searchlight Brain RDM Calculation: ');
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smoothing and upsampling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    FFX_permutation(Models, userOptions)
    fprintf('Stage 2 - Fixed Effects Analysis: ');
    toc
else
    % random effect test
    fprintf('Performing permutation tests over all subjects (random effect) to calculate p-values.');
    tic
    RFX_permutation(Models,userOptions)
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
    MEGFindCluster_source(Models, range, userOptions);
end
fprintf('Stage 3 - Spatiotemporal Clustering: ' );
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute cluster level p-values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
get_cluster_p_value(Models,userOptions);
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
    deleteDir(userOptions, Models);
end

%%%%%%%%%%%%%%%%%%%%
%% Sending an email %%
%%%%%%%%%%%%%%%%%%%%
if userOptions.recieveEmail
    setupInternet();
    setupEmail(userOptions.mailto);
end

