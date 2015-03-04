%  Recipe_MEG_source
%
% Cai Wingfield 9-2010
% Updated: Isma Zulfiqar 11-2012

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

toolboxRoot = '/imaging/ls02/toolbox/devel/toolbox'; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code
userOptions = projectOptions();

Models = rsa.constructModelRDMs(userOptions);

userOptions = rsa.meg.setMetadata_MEG(Models, userOptions); 
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%
nSubjects = userOptions.nSubjects;
for subject =1:nSubjects
    thisSubject = userOptions.subjectNames{subject};
    fprintf(['Reading MEG source solutions for subject number ' num2str(subject) ' of ' num2str(nSubjects) ': ' thisSubject ':']);    
    sourceMeshes.(thisSubject) = rsa.meg.MEGDataPreparation_source(subject,userOptions.betaCorrespondence, userOptions);
end
indexMasks   = rsa.meg.MEGMaskPreparation_source(userOptions);
maskedMeshes = rsa.meg.MEGDataMasking_source(sourceMeshes, indexMasks, userOptions.betaCorrespondence, userOptions);

%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

RDMs  = rsa.constructRDMs(maskedMeshes, userOptions.betaCorrespondence, userOptions);
RDMs  = rsa.rdm/averageRDMs_subjectSession(RDMs, 'session');
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

% fixed effects analysis
rsa.stat.testSignificance({aRDMs}, {Models}, userOptions);

% random effects analysis
rsa.stat.testSignificance_RandomEffects(RDMs, Models, userOptions)
