% DEMO1_RSA_ROI_sim is a script.  It is a sample recipe file which will
% simulate some fMRI data and then run through an "RoI-based" RSA pipeline for
% the simulated data
%
% Cai Wingfield 5-2010, 6-2010, 7-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
clear;clc;
cd ..;
toolboxRoot = pwd; addpath(genpath(toolboxRoot));
cd Demos;
% Generate a userOptions structure and then clone it for the two streams of
% data
% in this pipeline. Change only the names.
mkdir DEMO2;
userOptions_common = projectOptions_demo();
% Generate a simulationOptions structure.
simulationOptions = simulationOptions_demo();

%%%%%%%%%%%%%%%%
%% Simulation %%
%%%%%%%%%%%%%%%%

% Generate the SPM files for each subject containing conditions clustered
% according to preferences in the simulationOptions.
promptOptions.checkFiles(1).address = fullfile(userOptions_common.rootPath, 'Details', [userOptions_common.analysisName, '_simulateDataFiles_Details.mat']);
promptOptions.checkFiles(2).address = fullfile(userOptions_common.rootPath, 'Details', [userOptions_common.analysisName, 'True_fMRIDataPreparation_Details.mat']);
promptOptions.checkFiles(3).address = fullfile(userOptions_common.rootPath, 'Details', [userOptions_common.analysisName, 'Noisy_fMRIDataPreparation_Details.mat']);

for fileCount = 1:numel(promptOptions.checkFiles)
    if exist(promptOptions.checkFiles(fileCount).address, 'file')
        promptFlag(fileCount) = true; % If any of the passed-in files already exist, this will get toggled to true
    end%if
end%for
if ~exist('promptFlag','var'), promptFlag = false;end
if prod(double(promptFlag))
    reply = input('You have done the simulations before.Do you want to run it again? Y/N [Y]: ', 's');
    if ismember(reply,{'n','N','NO','No','no'})
        userOptions_common.forcePromptReply = 'S';
    else
        userOptions_common.forcePromptReply = 'R';
    end  
end
userOptions_true = userOptions_common; userOptions_true.analysisName = [userOptions_true.analysisName 'True'];
userOptions_noisy = userOptions_common; userOptions_noisy.analysisName = [userOptions_noisy.analysisName 'Noisy'];
[betaCorrespondence_true,betaCorrespondence_noisy,fMRI] = rsa.sim.simulateDataFiles(userOptions_common, simulationOptions);

% Load in the 'true' fMRI data
fullBrainVols_true = rsa.fmri.fMRIDataPreparation(betaCorrespondence_true, userOptions_true);

% Load in the 'noisy' fMRI data
fullBrainVols_noisy = rsa.fmri.fMRIDataPreparation(betaCorrespondence_noisy, userOptions_noisy);

% Name the RoIs for both streams of data
RoIName = 'SimRoI';
responsePatterns_true.(['true' RoIName]) = fullBrainVols_true;
responsePatterns_noisy.(['noisy' RoIName]) = fullBrainVols_noisy;

%%%%%%%%%%
%% RDMs %%
%%%%%%%%%%

% Construct RDMs for the 'true' data. One RDM for each subject (sessions
% have
% not been simulated) and one for the average across subjects.
RDMs_true = rsa.constructRDMs(responsePatterns_true, betaCorrespondence_true, userOptions_true);
RDMs_true = rsa.rdm.averageRDMs_subjectSession(RDMs_true, 'session');
averageRDMs_true = rsa.rdm.averageRDMs_subjectSession(RDMs_true, 'subject');

% Do the same for the 'noisy' data.
RDMs_noisy = rsa.constructRDMs(responsePatterns_noisy, betaCorrespondence_noisy, userOptions_noisy);
RDMs_noisy = rsa.rdm.averageRDMs_subjectSession(RDMs_noisy, 'session');
averageRDMs_noisy = rsa.rdm.averageRDMs_subjectSession(RDMs_noisy, 'subject');

% Prepare the model RDMs.
RDMs_model = rsa.constructModelRDMs(modelRDMs_demo2, userOptions_common);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Display the three sets of RDMs: true, noisy and model
rsa.figureRDMs(rsa.rdm.concatenateRDMs(RDMs_true, averageRDMs_true), userOptions_true, struct('fileName', 'noiselessRDMs', 'figureNumber', 1));
rsa.figureRDMs(rsa.rdm.concatenateRDMs(RDMs_noisy, averageRDMs_noisy), userOptions_noisy, struct('fileName', 'noisyRDMs', 'figureNumber', 2));
rsa.figureRDMs(RDMs_model, userOptions_common, struct('fileName', 'modelRDMs', 'figureNumber', 3));
% 
% Determine dendrograms for the clustering of the conditions for the two data
% streams
[blankConditionLabels{1:size(RDMs_model(1).RDM, 2)}] = deal(' ');
rsa.dendrogramConditions(averageRDMs_true, userOptions_true, struct('titleString', 'Dendrogram of conditions without simulated noise', 'useAlternativeConditionLabels', true, 'alternativeConditionLabels', {blankConditionLabels}, 'figureNumber', 4));
rsa.dendrogramConditions(averageRDMs_noisy, userOptions_noisy, struct('titleString', 'Dendrogram of conditions with simulated noise', 'useAlternativeConditionLabels', true, 'alternativeConditionLabels', {blankConditionLabels}, 'figureNumber', 5));
% 
% Display MDS plots for the condition sets for both streams of data
rsa.MDSConditions(averageRDMs_true, userOptions_true, struct('titleString', 'MDS of conditions without simulated noise', 'alternativeConditionLabels', {blankConditionLabels}, 'figureNumber', 6));
rsa.MDSConditions(averageRDMs_noisy, userOptions_noisy, struct('titleString', 'MDS of conditions with simulated noise', 'alternativeConditionLabels', {blankConditionLabels}, 'figureNumber', 7));
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second-order analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display a second-order simmilarity matrix for the models and the true and noisy simulated pattern RDMs
rsa.pairwiseCorrelateRDMs({averageRDMs_true, averageRDMs_noisy, RDMs_model}, userOptions_common, struct('figureNumber', 8));

% Plot all RDMs on a MDS plot to visualise pairwise distances.
rsa.MDSRDMs({averageRDMs_true, averageRDMs_noisy, RDMs_model}, userOptions_common, struct('titleString', 'MDS of noisy RDMs and models', 'figureNumber', 11));
 
for i=1:numel(RDMs_model)
    models{i}=RDMs_model(i);
end
models{end+1} = averageRDMs_true;

%% statistical inference:
% test the relatedness and compare the candidate RDMs

userOptions = userOptions_noisy;
userOptions.RDMcorrelationType='Kendall_taua';
userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
userOptions.RDMrelatednessThreshold = 0.05;
userOptions.figureIndex = [10 11];
userOptions.RDMrelatednessMultipleTesting = 'FDR';
userOptions.candRDMdifferencesTest = 'subjectRFXsignedRank';
userOptions.candRDMdifferencesThreshold = 0.05;
userOptions.candRDMdifferencesMultipleTesting = 'none';
stats_p_r=rsa.compareRefRDM2candRDMs(RDMs_noisy, models, userOptions);
