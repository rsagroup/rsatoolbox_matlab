%  Recipe_MEG_source
%
% Cai Wingfield 9-2010
% Updated: Isma Zulfiqar 11-2012

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = '/imaging/ls02/toolbox/devel/toolbox'; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code
userOptions = projectOptions();

Models = constructModelRDMs(userOptions);

userOptions = setMetadata_MEG(Models, userOptions); 
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%
nSubjects = userOptions.nSubjects;
for subject =1:nSubjects
    thisSubject = userOptions.subjectNames{subject};
    fprintf(['Reading MEG source solutions for subject number ' num2str(subject) ' of ' num2str(nSubjects) ': ' thisSubject ':']);    
    sourceMeshes.(thisSubject) = MEGDataPreparation_source(subject,userOptions.betaCorrespondence, userOptions);
end
indexMasks = MEGMaskPreparation_source(userOptions);
maskedMeshes = MEGDataMasking_source(sourceMeshes, indexMasks, userOptions.betaCorrespondence, userOptions);

%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

RDMs = constructRDMs(maskedMeshes, userOptions.betaCorrespondence, userOptions);
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

% fixed effects analysis
testSignificance({aRDMs}, {Models}, userOptions);

% random effects analysis
testSignificance_RandomEffects(RDMs, Models, userOptions)
