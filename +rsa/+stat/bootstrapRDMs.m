function [realRs, bootstrapEs, pairwisePs, bootstrapRs] = bootstrapRDMs(bootstrappableReferenceRDMs, candRDMs, userOptions)
% [realRs bootstrapEs pairwisePs bootstrapRs] = ...
%                        bootstrapRDMs(bootstrappableReferenceRDMs, ...
%                                                candRDMs, ...
%                                                userOptions ...
%                                                )
%
%        bootstrappableReferenceRDMs --- The RDMs to bootstrap.
%                bootstrappableReferenceRDMs should be a [nConditions nConditions
%                nSubjects]-sized matrix of stacked squareform RDMs.
%
%        candRDMs --- The RDMs to test against.
%                testRDM should be an [nConditions nConditions nCandRDMs]-sized
%                matrix where each leaf is one RDM to be tested.
%
%        userOptions --- The options struct.
%                userOptions.corrType
%                        A string descriptive of the distance measure to be
%                        used. Defaults to 'Spearman'.
%                userOptions.nBootstrap
%                        How many bootstrap resamplings shoule be performed?
%                        Defaults to 1000.
%                userOptions.resampleSubjects
%                        Boolean. If true, subjects will be bootstrap resampled.
%                        Defaults to false.
%                userOptions.resampleConditions
%                        Boolean. If true, conditions will be resampled.
%                        Defaults to true.
%
%        realRs
%                The true RDM correlations between the average of the
%                bootstrappableReferenceRDMs and the testRDM.
%
%        bootstrapStdE
%                The bootstrap standard error.
%
% Cai Wingfield 6-2010, 7-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

% Sort out defaults
userOptions = setIfUnset(userOptions, 'nBootstrap', 1000);
userOptions = setIfUnset(userOptions, 'resampleSubjects', false);
userOptions = setIfUnset(userOptions, 'resampleConditions', true);
userOptions = setIfUnset(userOptions, 'distanceMeasure', 'Spearman');

% Constants
nConditions = size(bootstrappableReferenceRDMs, 1);
nSubjects = size(bootstrappableReferenceRDMs, 3);

nCandRDMs = size(candRDMs, 3);

if ~(size(bootstrappableReferenceRDMs, 1) == size(candRDMs, 1))
    error('bootstrapRDMComparison:DifferentSizedRDMs', 'Two RDMs being compared are of different sizes. This is incompatible\nwith bootstrap methods!');
end%if

averageReferenceRDM = sum(bootstrappableReferenceRDMs, 3) ./ nSubjects;

% Decide what to say
if userOptions.resampleSubjects
    if userOptions.resampleConditions
        message = 'subjects and conditions';
    else
        message = 'subjects';
    end%if
else
    if userOptions.resampleConditions
        message = 'conditions';
    else
        message = 'nothing';
        warning('(!) You''ve gotta resample something, else the bar graph won''t mean anything!');
    end%if
end%if

fprintf(['Resampling ' message ' ' num2str(userOptions.nBootstrap) ' times']);

tic; %1

% Come up with the random samples (with replacement)
if userOptions.resampleSubjects
    resampledSubjectIs = ceil(nSubjects * rand(userOptions.nBootstrap, nSubjects));
else
    resampledSubjectIs = repmat(1:nSubjects, userOptions.nBootstrap, nSubjects);
end%if:resampleSubjects

if userOptions.resampleConditions
    resampledConditionIs = ceil(nConditions * rand(userOptions.nBootstrap, nConditions));
else
    resampledConditionIs = repmat(1:nConditions, userOptions.nBootstrap, 1);
end%if:resampleConditions

% Preallocation
realRs = nan(nCandRDMs, 1);
bootstrapRs = nan(nCandRDMs, userOptions.nBootstrap);
bootstrapEs = nan(nCandRDMs, 1);
pairwisePs = nan(nCandRDMs, nCandRDMs);
% replace the diagonals for each instance of the candidate RDMs with
% NaN entries
for subI = 1:size(bootstrappableReferenceRDMs,3)
    temp = bootstrappableReferenceRDMs(:,:,subI);
    temp(logical(eye(size(temp,1)))) = nan;
    bootstrappableReferenceRDMs(:,:,subI) = temp;
end

% Bootstrap
n = 0; k=0;

fprintf('\n');

% Need to create one candidate RDM for each reference because we don't want
% to correlate averaged RDMs.
candRDMs3rd = cell(nCandRDMs, 1);
for candRDMI = 1:nCandRDMs
    candRDMs3rd{candRDMI} = repmat(candRDMs(:,:,candRDMI), ...
        [1, 1, nSubjects]);
end

for candRDMI = 1:nCandRDMs
    for b = 1:userOptions.nBootstrap
        n = n + 1;
        localReferenceRDMs = bootstrappableReferenceRDMs(resampledConditionIs(b,:),resampledConditionIs(b,:),resampledSubjectIs(b,:));
        localTestRDM = candRDMs3rd{candRDMI}(resampledConditionIs(b,:), resampledConditionIs(b,:), resampledSubjectIs(b,:));

        if isequal(userOptions.RDMcorrelationType,'Kendall_taua')
            %             tic
            bootstrapRs(candRDMI, b)=mean(diag(rankCorr_Kendall_taua(squeeze(vectorizeRDMs(localReferenceRDMs)),squeeze(vectorizeRDMs(localTestRDM)))));
            %             toc
        elseif isequal(userOptions.RDMcorrelationType,'raeSpearman')
            %             tic
            %             bootstrapRs(candRDMI, b)=raeSpearmanCorr(vectorizeRDMs(averageBootstrappedRDM),vectorizeRDMs(localTestRDM));
            %             toc
            %             tic
            bootstrapRs(candRDMI, b)=mean(diag(raeSpearmanCorr(vectorizeRDMs(localReferenceRDMs),vectorizeRDMs(localTestRDM))));
            %             toc
        else
            %             tic
            bootstrapRs(candRDMI, b) = mean(diag(corr(squeeze(vectorizeRDMs(localReferenceRDMs)), squeeze(vectorizeRDMs(localTestRDM)), 'type',userOptions.distanceMeasure,'rows','pairwise')));
            %             toc
        end

        if mod(n,floor(userOptions.nBootstrap*nCandRDMs/100))==0
            fprintf('%d%% ',floor(100*n/(userOptions.nBootstrap*nCandRDMs)));
            if mod(n,floor(userOptions.nBootstrap*nCandRDMs/10))==0, fprintf('\n'); end;
        end
    end%for:b
end%for:candRDMI
fprintf('\n');

bootstrapEs = std(bootstrapRs, 0, 2);
k=0;
for candRDMI = 1:nCandRDMs
%     if isequal(userOptions.RDMcorrelationType,'Kendall_taua')
%         realRs(candRDMI)=rankCorr_Kendall_taua(vectorizeRDMs(averageReferenceRDM)',vectorizeRDMs(candRDMs(:,:,candRDMI))');
%     else
%         realRs(candRDMI) = corr(vectorizeRDMs(averageReferenceRDM)',vectorizeRDMs(candRDMs(:,:,candRDMI))','type',userOptions.distanceMeasure,'rows','pairwise');
%     end
    for candRDMJ = 1:nCandRDMs
        if candRDMI == candRDMJ
            pairwisePs(candRDMI, candRDMJ) = nan;
        else
            ijDifferences = bootstrapRs(candRDMI, :) - bootstrapRs(candRDMJ, :);
            pairwisePs(candRDMI, candRDMJ) = 2*min([numel(find(ijDifferences < 0)),numel(find(ijDifferences > 0))] / userOptions.nBootstrap);
        end%if:diagonal
    end%for:candRDMJ
end%for:candRDMI

t = toc;%1

fprintf([': [' num2str(ceil(t)) 's]\n']);

end%function
