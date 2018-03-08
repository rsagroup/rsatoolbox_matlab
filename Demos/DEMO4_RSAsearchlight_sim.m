%% DEMO4_RSAsearchlight_sim
% simulates fMRI data for a number of subjects and runs searchlight
% analysis using RSA, computes the similarity maps for each subject
% and does group level inference.

% demo of the searchlight analysis
%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
clear;clc
returnHere = pwd; % We'll come back here later
cd ..;
toolboxRoot = pwd; addpath(genpath(toolboxRoot));
cd Demos;
mkdir('DEMO4');
% Generate a userOptions structure
userOptions = projectOptions_demo();
userOptions.rootPath = [pwd,filesep,'DEMO4'];
userOptions.analysisName = 'DEMO4';

% Generate a simulationOptions structure.
simulationOptions = simulationOptions_demo_SL();

searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;
searchlightOptions.nConditions = 40;
load([returnHere,filesep,'sampleMask_org.mat']);
load([returnHere,filesep,'anatomy.mat']);% load the resliced structural image

models = rsa.constructModelRDMs(modelRDMs_SL_sim, userOptions);

nCond = searchlightOptions.nConditions;
Nsubjects = 20;

warpFlags.interp = 1;
warpFlags.wrap = [0 0 0];
userOptions.voxelSize = [3 3 3];
warpFlags.vox = userOptions.voxelSize; % [3 3 3.75]
warpFlags.bb = [-78 -112 -50; 78 76 85];
warpFlags.preserve = 0;
mapsFilename = [userOptions.analysisName, '_fMRISearchlight_Maps.mat'];
RDMsFilename = [userOptions.analysisName, '_fMRISearchlight_RDMs.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRISearchlight_Details.mat'];

%% simulate the data and compute the correlation maps per subject

for subI = 1:Nsubjects
    subject = ['subject',num2str(subI)];
    maskName = 'mask';
    fprintf(['simulating fullBrain volumes for subject %d \n'],subI)
    
    [B_true,Mask,Y_true, fMRI_sub] = rsa.sim.simulateClusteredfMRIData_fullBrain(simulationOptions);
    B_noisy = fMRI_sub.B;
    singleSubjectVols = B_noisy';
    userOptions.searchlightRadius = 9;mask = m;
    fprintf(['computing correlation maps for subject %d \n'],subI)
    [rs, ps, ns, searchlightRDMs.(subject)] = rsa.fmri.searchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions);
    rsa.util.gotoDir(userOptions.rootPath, 'Maps');
    save(['rs_',subject,'.mat'],'rs');
    clear rs searchlightRDMs;
    cd(returnHere);
end
%% display the design matrix, model RDMs and simulated RDMs and the SL
rsa.fig.selectPlot(1);
subplot(321);imagesc(fMRI_sub.X);axis square
ylabel('\bfscans');xlabel('\bfregressors')
title('\bfdesign matrix')

subplot(322);plot(fMRI_sub.X(:,12:14))
xlabel('scans');title('\bfregressors for 3 example conditions')

subplot(323);
image(rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(squareform(pdist(fMRI_sub.groundTruth)),1)),'CDataMapping','scaled','AlphaData',~isnan(squareform(pdist(fMRI_sub.groundTruth))));
axis square off
title('\bfsimulated ground truth RDM')

subplot(324);
image(rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(models(1).RDM,1)),'CDataMapping','scaled','AlphaData',~isnan(models(1).RDM));
axis square off
colormap(rsa.fig.RDMcolormap)
title('\bftested model RDM')


relRoi = rsa.fmri.sphericalRelativeRoi(userOptions.searchlightRadius,userOptions.voxelSize);
nVox_searchlight = size(relRoi,1);
rsa.fig.showVoxObj(relRoi+repmat(simulationOptions.effectCen,[nVox_searchlight,1]),1,[3 2 5]);
title(['\bf searchlight with ',num2str(nVox_searchlight),' voxels'])
rsa.fig.handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'SLsimulationSettings'],userOptions);

%% load the previously computed rMaps and concatenate across subjects
% prepare the rMaps:
for subjectI = 1:Nsubjects
    load([userOptions.rootPath,filesep,'Maps',filesep,'rs_subject',num2str(subjectI),'.mat']);
    rMaps{subjectI} = rs;
    fprintf(['loading the correlation maps for subject %d \n'],subjectI);
end
% concatenate across subjects
for modelI = 1:numel(models)
    for subI = 1:Nsubjects
        thisRs = rMaps{subI};
        thisModelSims(:,:,:,subI) = thisRs(:,:,:,modelI);
    end
    % obtain a pMaps from applying a 1-sided signrank test and also t-test to
    % the model similarities:
    for x=1:size(thisModelSims,1)
        for y=1:size(thisModelSims,2)
            for z=1:size(thisModelSims,3)
                if mask(x,y,z) == 1
                    [h p1(x,y,z)] = ttest(squeeze(thisModelSims(x,y,z,:)),0,0.05,'right');
                    [p2(x,y,z)] = rsa.util.signrank_onesided(squeeze(thisModelSims(x,y,z,:)));
                else
                    p1(x,y,z) = NaN;
                    p2(x,y,z) = NaN;
                end
            end
        end
        disp(x);
    end
    % apply FDR correction
    pThrsh_t  = rsa.stat.FDRthreshold(p1,0.05,mask);
    pThrsh_sr = rsa.stat.FDRthreshold(p2,0.05,mask);
    p_bnf = 0.05/sum(mask(:));
    % mark the suprathreshold voxels in yellow
    supraThreshMarked_t = zeros(size(p1));
    supraThreshMarked_t(p1 <= pThrsh_t) = 1;
    supraThreshMarked_sr = zeros(size(p2));
    supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;
    
    % display the location where the effect was inserted (in green):
    brainVol = rsa.gmri.addRoiToVol(map2vol(anatVol),mask2roi(mask),[1 0 0],2);
    brainVol_effectLoc = rsa.gmri.addBinaryMapToVol(brainVol,Mask.*mask,[0 1 0]);
    rsa.fig.showVol(brainVol_effectLoc,'simulated effect [green]',2);
    rsa.fig.handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'results_DEMO4_simulatedEffectRegion'],userOptions);
    
    % display the FDR-thresholded maps on a sample anatomy (signed rank test) :
    brainVol = rsa.fmri.addRoiToVol(map2vol(anatVol),mask2roi(mask),[1 0 0],2);
    brainVol_sr = rsa.fmri.addBinaryMapToVol(brainVol,supraThreshMarked_sr.*mask,[1 1 0]);
    rsa.fig.showVol(brainVol_sr,'signrank, E(FDR) < .05',3)
    rsa.fig.handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'results_DEMO4_signRank'],userOptions);
    
    % display the FDR-thresholded maps on a sample anatomy (t-test) :
    brainVol = rsa.fmri.addRoiToVol(map2vol(anatVol),mask2roi(mask),[1 0 0],2);
    brainVol_t = rsa.fmri.addBinaryMapToVol(brainVol,supraThreshMarked_t.*mask,[1 1 0]);
    rsa.fig.showVol(brainVol_t,'t-test, E(FDR) < .05',4)
    rsa.fig.handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'results_DEMO2_tTest'],userOptions);
end

cd(returnHere);
