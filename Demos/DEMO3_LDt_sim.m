%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
clear;clc
cd ..
toolboxRoot = pwd; addpath(genpath(toolboxRoot));
cd Demos
mkdir('DEMO3')
% Generate a userOptions structure and then clone it for the two streams of data
% in this pipeline. Change only the names.
userOptions = projectOptions_demo();
userOptions.rootPath = [pwd,filesep,'DEMO3'];
userOptions.analysisName = 'DEMO3';

% Generate a simulationOptions structure.
simulationOptions = simulationOptions_demo_LDt();
%% simulate data for two sessions per subject and compute the
%% subject-specific LDAtRDMs:
for subI = 1:20%numel(userOptions_common.subjectNames)
    fprintf(['simulating data and computing LD-t values for subject %d \n'],subI)
    [B_true,Y_true,fMRI_a,fMRI_b] = rsa.sim.simulateClusteredfMRIData(simulationOptions);
    [RDM_fdtFolded_ltv, cv2RDM_fdt_sq] = rsa.stat.fisherDiscrTRDM(fMRI_a.X,fMRI_a.Y,fMRI_b.X,fMRI_b.Y);
    RDM_lda = squareform(RDM_fdtFolded_ltv);% diagonals will contain zeros
    RDMs(subI).RDM = RDM_lda;
    RDMs(subI).name = ['LDAtRDM | subject',num2str(subI)];
    RDMs(subI).color = [1 0 0]; 
end
%% compute the subject-averaged LDAtRDM
averageRDMs_LDt = rsa.rdm.averageRDMs_subjectSession(RDMs, 'subject');
% %% display RDMs for all subjects
% showRDMs(RDMs,1)
%% display the average LDAtRDM
averageRDMs_LDA.name = 'subject-average LD-tRDM';
rsa.fig.showRDMs(averageRDMs_LDt,2)
filespec = 'demo_LDAtRDM_simulation_1';
rsa.fig.handleCurrentFigure([userOptions.rootPath,filesep,'demo_LDAtRDM_simulation_groupAverageRDM'],userOptions);
    
%% random effects analysis
rdms = rsa.rdm.unwrapRDMs(RDMs);% nCond x nCond x nSubj
nCond = size(rdms,1);
p_t = ones(nCond,nCond);
p_sr = ones(nCond,nCond);
for condI = 1:nCond
    for condJ = condI+1:nCond
        [h p_t(condI,condJ)] = ttest(squeeze(rdms(condI,condJ,:)),0,0.05,'right');
        [p_sr(condI,condJ)] = rsa.stat.signrank_onesided(squeeze(rdms(condI,condJ,:)));
    end
end
%% compute the thresholds and display
thresh_uncorr = 0.05;
nTests = nCond*(nCond-1)/2;
thresh_fdr_t  = rsa.stat.FDRthreshold(p_t,thresh_uncorr);
thresh_fdr_sr = rsa.stat.FDRthreshold(p_sr,thresh_uncorr);
thresh_bnf = thresh_uncorr/nTests;
rsa.fig.selectPlot(5);
subplot(131);rsa.fig.image_thr(p_sr,thresh_uncorr);
axis square off;
title('\bfp < 0.05 (uncorr.)');

subplot(132);rsa.fig.image_thr(p_sr,thresh_fdr_t);
axis square off;
title('\bfp < 0.05 (FDR)');

subplot(133);rsa.fig.image_thr(p_sr,thresh_bnf);
axis square off;
title('\bfp < 0.05 (Bonferroni)');

filespec = 'demo_LDAtRDM_simulation_2';
rsa.fig.addHeading('random effect analysis, subjects as random effects')
rsa.fig.handleCurrentFigure([userOptions.rootPath,filesep,'demo_LDAtRDM_simulation_subjectRFX'],userOptions);
