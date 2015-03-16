function [ceiling_upperBound, ceiling_lowerBound, bestFitRDM]=ceilingAvgRDMcorr(refRDMestimates,RDMcorrelationType,monitor)

% Given a set of reference RDM estimates (e.g. from multiple subjects) in
% argument refRDMestimates, this function estimates upper and lower bounds
% on the ceiling, i.e. the highest average RDM correlation (across the RDM
% estimates) that the true model's RDM prediction achieves given the
% variability of the estimates.
%
% The upper bound on the ceiling is estimated by finding the hypothetical
% model RDM that maximises the average correlation to the reference RDM
% estimates (using the correlation coefficient specified by argument
% RDMcorrelationType).
%
% For Pearson and Spearman correlation, we have closed-form solutions. For
% Pearson correlation, the exact upper bound is the average correlation
% across subjects to the mean RDM computed after z-transforming the
% dissimilarities in each subject's RDM. For Spearman correlation, the
% exact upper bound is the average correlation across subjects to the mean
% RDM computed after rank-transforming the dissimilarities in each
% subject's RDM.
%
% For Kandall's tau a, there is a closed-form solution for an upper bound
% on the ceiling utilising the fact that the tau a correlation distance is
% proportional to the squared Euclidean distance in pair-relation space.
% Unfortunately, this upper bound is not tight enough to be useful. We
% therefore initially estimate the upper bound as the average Kendall tau a
% between the mean RDM computed after rank-transforming each RDM, and each
% subject's RDM. This provides a slight underestimate of the upper bound.
% We therefore attempt to iteratively optimise the estimate by searching an
% RDM that yields a higher average tau a correlation with the
% single-subject RDMs. In practice, we have observed only minimal
% improvements. When the RDM space is very high-dimensional, search is
% difficult and tau a takes long to compute, we therefore terminate the
% search after historyLength many iterations and return the best estimate
% of the upper bound (which slightly underestimates the true upper bound).
%
% We estimate the lower bound of the ceiling using a leave-one-subject-out
% crossvalidation procedure.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% preparations
if ~exist('RDMcorrelationType','var'), RDMcorrelationType='Kendall_taua'; end
if ~exist('monitor','var'), monitor=false; end
disp('Estimating the ceiling for the average RDM correlation...');

refRDMestimates=vectorizeRDMs(unwrapRDMs(refRDMestimates));
[ignore nDissims nSubjects]=size(refRDMestimates);


%% compute the upper bound for the ceiling
switch  RDMcorrelationType
    case 'Pearson'
        % Pearson correlation distance = const * Euclidean dist. squared in
        % z-transformed RDM space. Thus, z-transform RDMs to make the mean
        % RDM minimise the average correlation distance to the
        % single-subject RDMs.
        refRDMestimates=refRDMestimates-repmat(mean(refRDMestimates,2),[1 nDissims 1]);
        refRDMestimates=refRDMestimates./repmat(std(refRDMestimates,[],2),[1 nDissims 1]);
        bestFitRDM=mean(refRDMestimates,3);
        
    case 'Spearman'
        % Spearman correlation distance = const * Euclidean dist squares in
        % rank-transformed RDM space. Thus, rank-transform RDMs to make the
        % mean RDM minimise the average correlation distance to the
        % single-subject RDMs.
        refRDMestimates=reshape(tiedrank(reshape(refRDMestimates,[nDissims nSubjects])),[1 nDissims nSubjects]);
        bestFitRDM=mean(refRDMestimates,3);
        
    case 'Kendall_taua'
        % Closed-form hard upper bound via Euclidean distance in pair-relation space
        %         nDissimPairs=(nDissims^2-nDissims)/2;
        %         meanOfPairRelations=zeros(1,nDissimPairs);
        %         for subjectI=1:nSubjects
        %             cSubj_RDM=single(refRDMestimates(1,:,subjectI));
        %             pairRelations=sign(repmat(cSubj_RDM,[nDissims 1])-repmat(cSubj_RDM',[1 nDissims]));
        %             meanOfPairRelations=meanOfPairRelations+squareform(pairRelations);
        %         end
        %         meanOfPairRelations=meanOfPairRelations/nSubjects;
        %
        %         squaredEuclideanDistInPairRelSpace=nan(nSubjects,1);
        %         for subjectI=1:nSubjects % recompute pair relations in another loop to save memory, because these can be huge.
        %             cSubj_RDM=single(refRDMestimates(1,:,subjectI));
        %             pairRelations=sign(repmat(cSubj_RDM,[nDissims 1])-repmat(cSubj_RDM',[1 nDissims]));
        %             squaredEuclideanDistInPairRelSpace(subjectI)=sum( (meanOfPairRelations-squareform(pairRelations)).^2 );
        %         end
        %         tauas=1-squaredEuclideanDistInPairRelSpace/(2*nDissimPairs);
        %         ceiling_upperBound=mean(tauas);
        %         % This is a closed-form hard upper bound.
        %         % Unfortunately, it is not tight at all.
        %         bestFitRDM=nan;
        
        % No closed-form solution providing a tight upper bound for the
        % ceiling (to our knowledge), so initialise with the mean of the
        % rank-transformed RDMs and optimise iteratively. 
        refRDMestimates=reshape(tiedrank(reshape(refRDMestimates,[nDissims nSubjects])),[1 nDissims nSubjects]);
        bestFitRDM=mean(refRDMestimates,3);
end    


%% estimate lower bound on the ceiling
for subjectI = 1:nSubjects
    cSubjectRDM=refRDMestimates(:,:,subjectI);
    
    refRDMs_LOO = refRDMestimates;
    refRDMs_LOO(:,:,subjectI) = [];
    avgRefRDM_LOO = nanmean(refRDMs_LOO,3);

    RDMcorrs_LOO(subjectI)=correlation(cSubjectRDM,avgRefRDM_LOO,RDMcorrelationType);
end
ceiling_lowerBound=mean(RDMcorrs_LOO);


%% estimate the upper bound on the ceiling
for subjectI=1:nSubjects
    meanRDM_corrs(subjectI)=correlation(bestFitRDM,refRDMestimates(1,:,subjectI),RDMcorrelationType);
end
meanRDM_avgCorr=mean(meanRDM_corrs);
ceiling_upperBound=meanRDM_avgCorr;

if isequal(RDMcorrelationType,'Pearson') || isequal(RDMcorrelationType,'Spearman')
    return; % done via the appropriate closed-form solution
end
    
% use only half the subjects to test convergence speed
% bestFitRDM=mean(refRDMestimates(1,:,logical( round(rand(nSubjects,1)) ) ),3);


%% find best-fit RDM by gradient descent to estimate upper bound for ceiling
disp('Iterative optimisation of the upper bound...');
nExploredDirections=2;
historyLength=20;

stepSize=1/nSubjects/10;
iteration=1;
bestFitRDMcorrHistory=nan(historyLength,1);
proportionOfSuccessfulSamplesHistory=nan(historyLength,1);
proportionOfStationarySamplesHistory=nan(historyLength,1);
stepSizeHistory=nan(historyLength,1);

convergenceZoneLength=5;
convergenceStdThreshold=0.0001;
converged=false;

subjectSimplexExploration=false;
cSubjectWeights=ones(nSubjects,1)/nSubjects; % initialise to the subject mean

% if matlabpool('size')==0
%     nLabs=12;
%     while true
%         try matlabpool(nLabs); break;
%         catch ME, nLabs=nLabs-2; end
%     end
% end

while true
    
    if subjectSimplexExploration
        % explore the simplex of the single-subject reference RDM estimates
        subjectWeights=repmat(cSubjectWeights,[1 nExploredDirections])+stepSize*randn(nSubjects,nExploredDirections);
        subjectWeights=subjectWeights-repmat(min(subjectWeights,[],1),[nSubjects 1]); % ensure weights>=0
        subjectWeights=subjectWeights./repmat(sum(subjectWeights,1),[nSubjects 1]); % ensure weights sum to 1
        candBestFitRDMs=reshape(squeeze(refRDMestimates)*subjectWeights,[1 nDissims nExploredDirections]); % 1 by nDissims by nExploredDirections
    else
        % RANDOM directions
        candBestFitRDMs=repmat(bestFitRDM,[1 1 nExploredDirections])+stepSize*randn(1,nDissims,nExploredDirections);
        %candBestFitRDMs(candBestFitRDMs<0)=0; % prevent exploration of negative dissimilarities
    end
    
    candBestFitRDMs=cat(3,candBestFitRDMs,bestFitRDM); % the last one is the current best-fit RDM
    RDMcorrs=nan(nExploredDirections,nSubjects);
    
    for subjectI=1:nSubjects
        %parfor directionI=1:nExploredDirections+1
        for directionI=1:nExploredDirections+1
            RDMcorrs(directionI,subjectI)=correlation(candBestFitRDMs(1,:,directionI),refRDMestimates(1,:,subjectI),RDMcorrelationType);
        end % direction I
        fprintf([num2str(round(((iteration-1)*nSubjects+subjectI)/(historyLength*nSubjects)*100)),' %%']);
    end % subjectI
    cBestFit_avgRDMcorr=mean(RDMcorrs(end,:),2);
    explDir_avgRDMcorrs=mean(RDMcorrs(1:end-1,:),2);
    bestFitRDMcorrHistory(iteration)=cBestFit_avgRDMcorr;

    % test for convergence
    if iteration>=convergenceZoneLength && std(bestFitRDMcorrHistory(iteration-convergenceZoneLength+1:iteration))<convergenceStdThreshold, converged=true; end

    % visualise
    if monitor && (mod(iteration,2)==0 || converged || iteration>=historyLength)
        h=figure(400); set(h,'Color','w'); clf;
        
        subplot(3,1,1); hold on;
        plot([1 iteration],[meanRDM_avgCorr meanRDM_avgCorr],'-k','LineWidth',5);
        %plot([1 iteration],[rankMeanRDM_avgCorr rankMeanRDM_avgCorr],':k','LineWidth',4);
        plot(bestFitRDMcorrHistory(1:iteration),'r','LineWidth',2); xlabel('time step'); ylabel('current best-fit RDM''s average corr. to ref. RDM estimates','LineWidth',2);
        legend('mean-RDM avg. corr.','optimised-RDM avg. corr.','Location','SouthEast');
        title(['current best-fit RDM''s average ',RDMcorrelationType,' correlation to the RDM estimates: ',num2str(cBestFit_avgRDMcorr)]);
        
        subplot(3,2,3); plot(stepSizeHistory(1:iteration),'LineWidth',3); xlabel('iteration'); ylabel('step size');
        subplot(3,2,4); hold on;
        plot(proportionOfSuccessfulSamplesHistory(1:iteration),'r','LineWidth',4); xlabel('iteration');
        plot(proportionOfStationarySamplesHistory(1:iteration),'b','LineWidth',2); 
        xlabel('iteration'); ylabel('proportions of samples');
        legend('successful','stationary','Location','NorthEast');
        
        showRDMs(bestFitRDM,[400 3 2 5]);
        subplot(3,2,5); title('current best-fit RDM');
        drawnow;
    end
    
    % break if converged or time out or subject pressed 'p' for 'proceed'
    if converged, 
        fprintf('(of max # iterations) ')
        break; 
    end
    
    if iteration>=historyLength
        warning('ceilingAvgRDMcorr: ceiling estimation procedure did not converge.');
        break;
    end
   
    % check if any of the explored direction led to an improvement
    avgRDMcorrDiffs=explDir_avgRDMcorrs-cBestFit_avgRDMcorr;
    
    [mxImprovement mxImprovementI]=max(avgRDMcorrDiffs);
    if mxImprovement>0
        bestFitRDM=candBestFitRDMs(1,:,mxImprovementI);
        if subjectSimplexExploration
            cSubjectWeights=subjectWeights(:,mxImprovementI);
        end
        disp('improved.');
    end
    
    proportionOfSuccessfulSamplesHistory(iteration)=sum(avgRDMcorrDiffs>0)/nExploredDirections;
    proportionOfStationarySamplesHistory(iteration)=sum(avgRDMcorrDiffs==0)/nExploredDirections;
    stepSizeHistory(iteration)=stepSize;
    
    if 0.3<sum(avgRDMcorrDiffs>0)/nExploredDirections || 0.1<sum(avgRDMcorrDiffs==0)/nExploredDirections
        % more than 30% of our samples improve the fit 
        % or more than 10% of our samples are on a local plateau,
        % either way: take bigger steps
        stepSize=stepSize*(2^(1/2)); % move faster
    elseif sum(avgRDMcorrDiffs>0)/nExploredDirections<0.1
        % less than 10% of the explored directions improve the fit
        % (and we're not on a plateau). might be close to the peak, where
        % every way is down, so look more locally.
        stepSize=stepSize/(2^(1/2)); % move slower
    end
    
    iteration=iteration+1;
end % while true

% define ceiling upper bound
ceiling_upperBound=cBestFit_avgRDMcorr;


end%function


%% compute correlations of all types
function r=correlation(a,b,RDMcorrelationType)
	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
	
    switch RDMcorrelationType
        case 'negEuclidean'
            r=-sqrt(sum((b(:)-a(:)).^2));
        case 'Kendall_taua'
            r=rankCorr_Kendall_taua(a(:),b(:));
        otherwise
            r=corr(a(:),b(:),'type',RDMcorrelationType);
    end

end%function