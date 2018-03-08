function stats_p_r=compareRefRDM2candRDMs(refRDM, candRDMs, userOptions)

% GENERAL PURPOSE
% This function compares a reference RDM to multiple candidate RDMs. For
% example, the reference RDM could be a brain region's RDM and the
% candidate RDMs could be multiple computational models. Alternatively, the
% reference RDM could be a model RDM and the candidate RDMs could be
% multiple brain regions' RDMs. More generally, the candidate RDMs could
% include both model RDMs and RDMs from brain regions, and the reference
% RDM could be either a brain RDM or a model RDM. In all these cases, one
% reference RDM is compared to multiple candidates.
%
% TESTING RDM CORRELATIONS
% The function compares the reference RDM to each of the candidates, tests
% each candidate for significant RDM correlation (test dependent on input
% data and userOptions) and presents a bar graph of the RDM correlations in
% descending order from left to right (best candidates on the left) by
% default, or in the order in which the candidate RDMs are passed. In
% addition, pairwise comparisons between the candidate RDMs are performed
% to assess, for each pair of candidate RDMs, which of the two RDMs is more
% consistent with the reference RDM. A significant pairwise difference is
% indicated by a horizontal line above the corresponding two bars. Each bar
% comes with an error bar, which indicates the standard error, estimated by
% the same procedure as is used for the pairwise candidate comparisons
% (dependent on input data and userOptions, see below).
%
% Statistical inference on the correlation between the reference RDM and
% each candidate RDM is performed using a one-sided Wilcoxon signed-rank
% across subjects by default. When the number of subjects is insufficient
% (<12), or when requested in userOptions, the test is instead performed by
% condition-label randomisation. By default, the false-discovery rate is
% controlled for these tests across candidate models. When requested in
% userOptions, the familywise error rate is controlled instead.
%
% COMPARISONS BETWEEN CANDIDATE RDMs
% For the comparisons between candidate RDMs as well, the inference
% procedure depends on the input data provided and on userOptions. By
% default, a signed-rank test across repeated measurements of the RDMs
% (e.g. multiple subjects or sessions) is performed. Alternatively, these
% tests are performed by bootstrapping of the subjects and/or conditions
% set. Across the multiple pairwise comparisons, the function controls the
% familywise error rate (Bonferroni method) or the false discovery rate
% (Benjamini and Hochberg, 1995).
%
% CEILING UPPER AND LOWER BOUNDS
% If multiple instances of the reference RDM are passed (typically
% representing estimates for multiple subjects), a ceiling estimate is
% indicated by a gray transparent horizontal bar. The ceiling is the
% expected value, given the noise (i.e. the variability across subjects),
% of the average correlation of the true model's RDM with the
% single-subject RDMs. The upper and lower edges of the ceiling bar are
% upper- and lower-bound estimates for the unknown true ceiling. The upper
% bound is estimated by computing the correlation between each subject's
% RDM and the group-average RDM. The lower bound is estimated similarly,
% but using a leave-one-out approach, where each subject's RDM is
% correlated with the average RDM of the other subjects' RDMs. When Pearson
% correlation is chosen for comparing RDMs
% (userOptions.RDMcorrelationType), the RDMs are first z-transformed. When
% Spearman correlation is chosen, the RDMs are first rank-transformed. When
% Kendall's tau a is chosen, an iterative procedure is used. See Nili et
% al. (2013) for a full motivation for the ceiling estimate and for an
% explanation of why these are useful upper- and lower-bound estimates.
%
% USAGE
%       stats_p_r=compareRefRDM2candRDMs(refRDM, candRDMs[, userOptions])
%
% ARGUMENTS
% refRDM
%       The reference RDM, which can be a square RDM or lower-triangular-
%       vector RDM, or a wrapped RDM (structure specifying a name and
%       colour for colour-coding of the RDM in figures). refRDM may also be
%       a set of independent estimates of the same RDM (square matrices or
%       lower-triangular vectors stacked along the third dimension or a
%       structured array of wrapped RDMs), e.g. an estimate for each of
%       multiple subjects or sessions, which are then used for
%       random-effects inference.
%
% candRDMs
%       A cell array with one cell for each candidate RDM. The candidate
%       RDMs can be square or lower-triangular-vector RDMs or wrapped RDMs.
%       Each candidate RDM may also contain multiple independent estimates
%       of the same RDM, e.g. an estimate for each of multiple subjects or
%       sessions. These can be used for random-effects inference if all
%       candidate RDMs have the same number of independent estimates,
%       greater than or equal to 12. However, if refRDM contains 12 or more
%       independent estimate of the reference RDMs, then random-effects
%       inference is based on these and multiple instances of any candidate
%       RDMs are replaced by their average. In case the dissimilarity for a
%       given pair of conditions is undefined (NaN) in any candidate RDM or
%       in the reference RDM, that pair is set to NaN in all RDMs and
%       treated as a missing value. This ensures that comparisons between
%       candidate RDMs are based on the same set of dissimilarities.
%
% userOptions.RDMcorrelationType
%       The correlation coefficient used to compare RDMs. This is
%       'Spearman' by default, because we prefer not to assume a linear
%       relationship between the distances (e.g. when a brain RDM from fMRI
%       is compared to an RDM predicted by a computational model).
%       Alternative definitions are 'Kendall_taua' (which is appropriate
%       whenever categorical models are tested) and 'Pearson'. The Pearson
%       correlation coefficient may be justified when RDMs from the same
%       origin (e.g. multiple computational models or multiple brain
%       regions measured with the same method) are compared. For more
%       details on the RDM correlation type, see Nili et al. (2013), in
%       particular, Figure S2 in the supplementary information.
%
% userOptions.RDMrelatednessTest
%       'subjectRFXsignedRank' (default): Test the relatedness of the reference RDM
%       to each candidate RDM by computing the correlation for each subject
%       and performing a one-sided Wilcoxon signed-rank test against the
%       null hypothesis of 0 correlation, so as to test for a positive
%       correlation. (Note that multiple independent measurements of the
%       reference or candidate RDMs could also come from repeated
%       measurements within one subject. We refer to the instances as
%       �subjects�, because subject random-effects inference is the most
%       common case.)
%
%       'randomisation': Test the relatedness of the reference RDM to each
%       candidate RDM by randomising the condition labels of the reference
%       RDM, so as to simulate the null distribution for the RDM
%       correlation with each candidate RDM. In case there are multiple
%       instances of the reference or candidate RDMs, these are first
%       averaged.
%
%       'conditionRFXbootstrap': Test the relatedness of the reference RDM
%       to each candidate RDM by bootstrapping the set of conditions
%       (typically: stimuli). For each bootstrap sample of the conditions
%       set, a new instance is generated for the reference RDM and for each
%       of the candidate RDMs. Because bootstrap resampling is resampling
%       with replacement, the same condition can appear multiple times in a
%       sample. This entails 0 entries (from the diagonal of the original
%       RDM) in off-diagonal positions of the RDM for a bootstrap sample.
%       These zeros are treated as missing values and excluded from the
%       dissimilarities, across which the RDM correlations are computed.
%       The p value for a one-sided test of the relatedness of each
%       candidate RDM to the reference RDM is computed as the proportion of
%       bootstrap samples with a zero or negative RDM correlation. This
%       test simulates the variability of the estimates across condition
%       samples and thus supports inference generalising to the population
%       of conditions (or stimuli) that the condition sample can be
%       considered a random sample of. Note that basic bootstrap tests are
%       known to be slightly optimistic.
%
%       'subjectConditionRFXbootstrap': Bootstrap resampling is
%       simultaneously performed across both subjects and conditions. This
%       simulates the greater amount of variability of the estimates
%       expected if the experiment were repeated with a different sample of
%       subjects and conditions. This more conservative test attempts to
%       support inference generalising across both subjects and stimuli (to
%       their respective populations). Again, the usual caveats for basic
%       bootstrap tests apply.
%
%       'none': Omit the test of RDM relatedness.
%
% userOptions.RDMrelatednessThreshold
%       The significance threshold (default: 0.05) for testing each
%       candidate RDM for relatedness to the reference RDM. Depending on
%       the choice of multiple testing correction (see next userOptions
%       field), this can be the expected false discovery rate, the
%       familywise error rate, or the uncorrected p threshold.
%
% userOptions.RDMrelatednessMultipleTesting
%       'FDR' (default): Control the false discovery rate across the
%       multiple tests (one for each candidate RDM). With this option,
%       userOptions.RDMrelatednessThreshold is interpreted to specify the
%       expected false discovery rate, i.e. the expected proportion of
%       candidate RDMs falsely declared significant among all candidate
%       RDMs declared significant.
%
%       'FWE': Control the familywise error rate. When the condition-label
%       randomisation procedure is selected to test RDM relatedness, then
%       randomisation is used to simulate the distribution of maximum RDM
%       correlations across all candidate RDMs. This method is more
%       powerful than Bonferroni correction when there are dependencies
%       among candidate RDMs. If another test is selected to test RDM
%       relatedness, the Bonferroni method is used. In either case,
%       userOptions.RDMrelatednessThreshold is interpreted as the
%       familywise error rate, i.e. the probability of getting any false
%       positives under the omnibus null hypothesis that all candidate RDMs
%       are unrelated to the reference RDM.
%
%       'none': Do not correct for multiple testing (not recommended). With
%       this setting, userOptions.RDMrelatednessThreshold is interpreted as
%       the uncorrected p threshold.
%
% userOptions.candRDMdifferencesTest
%       'subjectRFXsignedRank' (default, data permitting): For each pair of
%       candidate RDMs, perform a statistical comparison to determine which
%       candidate RDM better explains the reference RDM by using the
%       variability across subjects of the reference or candidate RDMs. The
%       test is a two-sided Wilcoxon signed-rank test of the null
%       hypothesis that the two RDM correlations (refRDM to each of the
%       candidate RDMs) are equal. This is the default test when multiple
%       instances of the reference RDM (typically corresponding to multiple
%       subjects) or a consistent number of multiple instances of each
%       candidate RDMs is provided and the number of multiple instances is
%       12 or greater. This test supports inference generalising to the
%       population of subjects (or repeated measurements) that the sample
%       can be considered a random sample of.
%
%       'subjectRFXbootstrap':  For each pair of candidate RDMs, perform a
%       two-sided statistical comparison to determine, which candidate RDM
%       better explains the reference RDM by bootstrapping the set of
%       subjects. For each bootstrap sample of the subjects set, the RDMs
%       are averaged across the bootstrap sample and the difference between
%       the two RDM correlations (refRDM to each of the candidate RDMs) is
%       computed. The p value is estimated as the proportion of bootstrap
%       samples further in the tails (symmetrically defined for a two-sided
%       test) than 0. This test simulates the variability of the estimates
%       across subject samples and thus supports inference generalising to
%       the population of subjects (or repeated measurements) that the
%       sample can be considered a random sample of. The usual caveats for
%       basic bootstrap tests apply.
%
%       'conditionRFXbootstrap': For each pair of candidate RDMs, perform a
%       two-sided statistical comparison to determine which candidate RDM
%       better explains the reference RDM by bootstrapping the set of
%       conditions (typically: stimuli). For each bootstrap sample of the
%       conditions set, a new instance is generated for the reference RDM
%       and for each of the candidate RDMs. Because bootstrap resampling is
%       is resampling with replacement, the same condition can appear
%       multiple times in a sample. This entails 0 entries (from the
%       diagonal of the original RDM) in off-diagonal positions of the RDM
%       for a bootstrap sample. These zeros are treated as missing values
%       and excluded from the dissimilarities, across which the RDM
%       correlations are computed. The p value for the two-sided test of
%       the difference for each pair of candidate RDMs is computed as for
%       the setting subjectRFXbootstrap (see above). This test simulates
%       the variability of the estimates across condition samples and thus
%       supports inference generalising to the population of conditions
%       (typically stimuli) that the condition sample can be considered a
%       random sample of. Again, the usual caveats for bootstrap tests
%       apply.
%
%       'subjectConditionRFXbootstrap': Bootstrap resampling is
%       simultaneously performed across both subjects and conditions. This
%       simulates the greater amount of variability of the estimates
%       expected if the experiment were repeated with a different sample of
%       subjects and conditions. This more conservative test attempts to
%       support inference generalising across both subjects and stimuli (to
%       their respective populations. However, the usual caveats for
%       bootstrap tests apply.
%
%       'none': Omit the pairwise tests of candidate RDMs comparing their
%       relative ability to explain the reference RDM.
%
% userOptions.candRDMdifferencesThreshold
%       The significance threshold for comparing each pair of candidate
%       RDMs in terms of their relatedness to the reference RDM. Depending
%       on the choice of multiple testing correction (see next userOptions
%       field), this can be the expected false discovery rate, the
%       familywise error rate, or the uncorrected p threshold.
%
% userOptions.candRDMdifferencesMultipleTesting
%       'FDR': Control the false discovery rate across the multiple tests
%       (one for each candidate RDM). With this option,
%       userOptions.candRDMdifferencesThreshold is interpreted to specify
%       the expected false discovery rate, i.e. the expected proportion of
%       pairs of candidate RDMs falsely declared significantly different
%       among all pairs of candidate RDMs declared significantly different
%       (in their relatedness to the reference RDM).
%
%       'FWE': Control the familywise error rate. With this option, the
%       Bonferroni method is used to ensure that the familywise error rate
%       is controlled. userOptions.candRDMdifferencesThreshold is
%       interpreted as the familywise error rate, i.e. the probability of
%       getting any false positives under the omnibus null hypothesis
%       that all pairs of candidate RDMs are equally related to the
%       reference RDM.
%
%       'none': Do not correct for multiple testing (not recommended).
%       userOptions.candRDMdifferencesThreshold is interpreted as the
%       uncorrected p threshold.
%
% userOptions.nRandomisations
%       The number of condition-label randomisations (default: 10,000) used
%       to simulate the null distribution that the reference RDM is
%       unrelated to each of the candidate RDMs.
%
% userOptions.nBootstrap
%       The number of bootstrap resamplings (default: 1,000) used in all
%       selected bootstrap procedures (relatedness test, candidate
%       comparison tests, error bars).
%
% userOptions.plotpValues
%       This option controls how the significance of the RDM relatedness
%       tests is indicated. If set to '*' (default), then an asterisk is
%       plotted on the bar for each candidate RDM that is significantly
%       related to the reference RDM. Significance depends on the test (see
%       above) and on userOptions.RDMrelatednessThreshold (default: 0.05)
%       and on userOptions.RDMrelatednessMultipleTesting (default: 'FDR').
%       Asterisks mark candidate RDMs that are significant at the specified
%       threshold and with the chosen method for accounting for multiple
%       testing. If set to '=', then the uncorrected p value is plotted
%       below the bar for each candidate RDM, and it is plotted in bold
%       type if it is significant by the criteria explained above.
%
% userOptions.barsOrderedByRDMCorr
%       This option controls the order of the displayed bars (default:
%       true). If set to true, bars corresponding to candidate RDMs are
%       displayed in descending order (from left to right) according to
%       their correlation to the reference RDM. Otherwise, bars are
%       displayed in the order in which the candidate RDMs are passed.
%
% userOptions.figureIndex (default: [1 2])
%       This option enables the user to specify the figure numbers for the
%       two created figures. The first figure contains the bargraph and the
%       second contains matrices indicating the significance of the
%       pairwise candidate RDM comparisons. The first panel shows the
%       uncorrected-p matrix. The second panel shows the thresholded
%       uncorrected-p matrix. The third panel shows the FDR-thresholded p
%       matrix. The fourth panel shows the Bonferroni-thresholded p matrix.
%
% userOptions.resultsPath  (default: pwd, i.e. current working directory)
%       This argument specifies the absolute path in which both figures are
%       to be saved.
%
% userOptions.saveFigurePDF
%       If set to true (default), the figures are saved in PDF format in
%       userOptions.resultsPath.
%
% userOptions.saveFigurePS
%       If true(default: false), the figures are saved in post-script
%       format in userOptions.resultsPath.
%
% userOptions.saveFigureFig
%       If true (default: false), the figures are saved in Matlab .fig
%       format in userOptions.resultsPath.
%
% userOptions.figure1filename
%       The filename for the bargraph figure, if chosen to be saved
%       (default: 'compareRefRDM2candRDMs_barGraph').
%
% userOptions.figure2filename
%       The filename for the p-value display figure, if chosen to be saved
%       (default: 'compareRefRDM2candRDMs_RDMcomparisonPvalues').
%
% RETURN VALUES
% stats_p_r
%       Structure containing numerical statistical results, including
%       effect sizes and p values.
%       stats_p_r.candRelatedness_r: average correlations to reference RDM
%       stats_p_r.candRelatedness_p: corresponding uncorrected p values
%       stats_p_r.SEs:               standard errors of average RDM
%                                    correlations
%       stats_p_r.candDifferences_r: matrix of bar-height differences
%                                    (i.e. average RDM-correlation
%                                    differences)
%       stats_p_r.candDifferences_p: matrix of p values for all pairwise
%                                    candidate comparisons
%       stats_p_r.orderedCandidateRDMnames: candidate RDM names in the
%                                    order in which the bars are displayed
%                                    (also the order used for the return
%                                    values)
%       stats_p_r.ceiling:           ceiling lower and upper bounds

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% set default options
userOptions = setIfUnset(userOptions, 'RDMcorrelationType', 'Spearman');
userOptions = setIfUnset(userOptions, 'RDMrelatednessTest', 'subjectRFXsignedRank');
userOptions = setIfUnset(userOptions, 'RDMrelatednessThreshold', 0.05);
userOptions = setIfUnset(userOptions, 'RDMrelatednessMultipleTesting', 'FDR');
userOptions = setIfUnset(userOptions, 'candRDMdifferencesTest','subjectRFXsignedRank');
userOptions = setIfUnset(userOptions, 'candRDMdifferencesThreshold',0.05);
userOptions = setIfUnset(userOptions, 'candRDMdifferencesMultipleTesting','FDR');
userOptions = setIfUnset(userOptions, 'plotpValues','*');
userOptions = setIfUnset(userOptions, 'nRandomisations',10000);
userOptions = setIfUnset(userOptions, 'nBootstrap',10000);
userOptions = setIfUnset(userOptions, 'barsOrderedByRDMCorr',true);
userOptions = setIfUnset(userOptions, 'figureIndex',[1 2]);
userOptions = setIfUnset(userOptions, 'resultsPath',pwd);
userOptions = setIfUnset(userOptions, 'saveFigurePDF', true);
userOptions = setIfUnset(userOptions, 'saveFigurePS', false);
userOptions = setIfUnset(userOptions, 'figure1filename', 'compareRefRDM2candRDMs_barGraph');
userOptions = setIfUnset(userOptions, 'figure2filename', 'compareRefRDM2candRDMs_RDMcomparisonPvalues');

%% inspect arguments and unwrap and average RDMs

refRDM_stack = stripNsquareRDMs(refRDM);
[nCond,nCond,nRefRDMinstances]=size(refRDM_stack);

if isstruct(refRDM)
   if isfield(refRDM,'name')
       refRDMName = refRDM.name;
   end
end
refRDMmask = ~any(isnan(refRDM_stack),3)|logical(eye(nCond)); % valid dissimilarities in the reference RDM
if any(~refRDMmask(:))
    nNaN_ref = numel(find(refRDMmask == 0));
    fprintf('Found %d NaNs in the reference RDM input. \n',nNaN_ref)
    fprintf('...excluding the entries from the analysis \n')
end
% collect the RDM labels in one cell array. also give names to unnamed
% RDMs.
unnamedI=1;
for candI = 1:numel(candRDMs)
    if isstruct(candRDMs{candI})
        try
            thisCand=candRDMs{candI};
            candidateNames{candI} = thisCand.name;
        catch
            candidateNames{candI} = ['unnamed RDM',num2str(unnamedI)];
            unnamedI = unnamedI + 1;
        end
    else
        candidateNames{candI} = ['unnamed RDM',num2str(unnamedI)];
        unnamedI = unnamedI + 1;
    end
    candRDM_stack{candI} = double(stripNsquareRDMs(candRDMs{candI}));
    meanCandRDMs(:,:,candI)= mean(candRDM_stack{candI},3);
    nsCandRDMinstances(candI)=size(candRDM_stack{candI},3);
end

candRDMmask = ~any(isnan(meanCandRDMs),3)|logical(eye(nCond));  % valid dissimilarities in all candidate RDMs

RDMmask = refRDMmask & candRDMmask; % valid dissimilarities in all RDMs
nValidDissimilarities=sum(RDMmask(:));

% find the indices for the non-NaN entries
validConds = ~any(~(RDMmask|logical(eye(nCond))),1);%~any(~(RDMmask|logical(eye(nCond))),1));

% reduce the condition set in both reference and candidate RDM stacks to
% locations where there are no NaNs only
refRDM_stack = refRDM_stack(validConds,validConds,:);
meanRefRDM = mean(refRDM_stack,3);

clear meanCandRDMs
for candI = 1:numel(candRDM_stack)
    temp=candRDM_stack{candI};
    candRDM_stack{candI}= temp(validConds,validConds,:);
    meanCandRDMs(:,:,candI) = mean(candRDM_stack{candI},3);
end


%% estimate the ceiling for the refRDM

nSubjects = size(refRDM_stack,3);
% for subj = 1:nSubjects
%     refrdms = refRDM_stack;
%     subjectRDM=refrdms(:,:,subj);
%     avgRefRDM = mean(refrdms,3);
%     refrdms(:,:,subj) = [];
%     avgRefRDM_LOO = mean(refrdms,3);
%     if isequal(userOptions.RDMcorrelationType,'Kendall_taua')
%         rs_ceilUpper(subj)=rankCorr_Kendall_taua(vectorizeRDMs(subjectRDM)',vectorizeRDMs(avgRefRDM)');
%         rs_ceilLower_LOO(subj)=rankCorr_Kendall_taua(vectorizeRDMs(subjectRDM)',vectorizeRDMs(avgRefRDM_LOO)');
%     elseif isequal(userOptions.RDMcorrelationType,'raeSpearman')
%         rs_ceilUpper(subj)=raeSpearmanCorr(vectorizeRDMs(subjectRDM)',vectorizeRDMs(avgRefRDM)');
%         rs_ceilLower_LOO(subj)=raeSpearmanCorr(vectorizeRDMs(subjectRDM)',vectorizeRDMs(avgRefRDM_LOO)');
%     else
%         rs_ceilUpper(subj)=corr(vectorizeRDMs(subjectRDM)',vectorizeRDMs(avgRefRDM)','type',userOptions.RDMcorrelationType,'rows','pairwise');
%         rs_ceilLower_LOO(subj)=corr(vectorizeRDMs(subjectRDM)',vectorizeRDMs(avgRefRDM_LOO)','type',userOptions.RDMcorrelationType,'rows','pairwise');
%     end
% end
% ceilingUpperBound=mean(rs_ceilUpper);
% ceilingLowerBound=mean(rs_ceilLower_LOO);
%
% stats_p_r.ceiling = [ceilingLowerBound,ceilingUpperBound];
%
[ceilingUpperBound, ceilingLowerBound, bestFitRDM]=ceilingAvgRDMcorr(refRDM_stack,userOptions.RDMcorrelationType,false);
stats_p_r.ceiling = [ceilingLowerBound,ceilingUpperBound];


%% display the average correlation bars
% correlate the average of the candRDMs with each instance of the refrDM
% cand2refSims: nRefRDMinstances x nCand
for candI = 1:numel(candRDMs)
    for subI = 1:nSubjects
        if isequal(userOptions.RDMcorrelationType,'Kendall_taua')
            cand2refSims(subI,candI)=rankCorr_Kendall_taua(vectorizeRDMs(meanCandRDMs(:,:,candI))',vectorizeRDMs(refRDM_stack(:,:,subI))');
        elseif isequal(userOptions.RDMcorrelationType,'raeSpearman')
            cand2refSims(subI,candI)=raeSpearmanCorr(vectorizeRDMs(meanCandRDMs(:,:,candI))',vectorizeRDMs(refRDM_stack(:,:,subI))');
        else
            cand2refSims(subI,candI)=corr(vectorizeRDMs(meanCandRDMs(:,:,candI))',vectorizeRDMs(refRDM_stack(:,:,subI))','type',userOptions.RDMcorrelationType,'rows','pairwise');
        end
    end
end
y = mean(cand2refSims,1);


[y_sorted,sortedIs]=sort(y,'descend');

for testi = 1:numel(candRDMs)
    for testj = 1:numel(candRDMs)
        stats_p_r.candDifferences_r(testi,testj) = y(testi)-y(testj);
    end
end

if ~userOptions.barsOrderedByRDMCorr
    sortedIs = 1:numel(candRDMs);
end


stats_p_r.candDifferences_r = stats_p_r.candDifferences_r(sortedIs,sortedIs);
candidateNames = candidateNames(sortedIs);
y_sorted = y(sortedIs);
stats_p_r.orderedCandidateRDMnames = candidateNames;
x = 1:numel(candRDMs);

h = figure(userOptions.figureIndex(1)); clf;
set(h, 'Color', 'w');
for barI = 1:numel(candRDMs)
    h=patch([barI-0.4 barI-0.4 barI+0.4 barI+0.4],[0 y_sorted(barI) y_sorted(barI) 0],[-0.01 -0.01 -0.01 -0.01],[68/225 131/225 149/225],'edgecolor','none');hold on;
    hold on
end

ylabel({'\bf RDM correlation ',['\rm[',deunderscore(userOptions.RDMcorrelationType),', averaged across ',num2str(nRefRDMinstances),' subjects]']});

stats_p_r.candRelatedness_r = cand2refSims(:,sortedIs);


%% decide inference procedures and test for RDM relatedness

% check if there are enough subjects for subject RFX tests
if nRefRDMinstances>=12
    subjectRFXacross='refRDMinstances';
    fprintf('Found %d instances of the reference RDM.',nRefRDMinstances);
    fprintf('Averaging all instances of each candidate RDM.\n');
else
    fprintf('Found less than 12 instances of the reference RDM.\n');

    if std(nsCandRDMinstances)==0 && nsCandRDMinstances(1)>=12
        fprintf('Found %d instances of all candidate RDMs. Using these for random-effects inference.\n',nsCandRDMinstances(1));
        fprintf('Averaging all instances of the reference RDM.\n');
        subjectRFXacross='candRDMinstances';
    else
        fprintf('Found less than 12 or inconsistent numbers of candidate-RDM instances. Cannot do subject random-effects inference.\n');
        subjectRFXacross='none';
    end
end

% check if there are enough conditions for condition RFX tests
if nCond >= 20
    conditionRFX=1;
else
    conditionRFX=0;
end

% decide RDM-relatedness test
if isequal(userOptions.RDMrelatednessTest,'subjectRFXsignedRank') || isequal(userOptions.RDMrelatednessTest,'subjectRFXbootstrap') || isequal(userOptions.RDMrelatednessTest,'subjectConditionRFXbootstrap')
    % user requested subject random-effects inference (by signed-rank or bootstrap test)

    if isequal(subjectRFXacross,'refRDMinstances');
        fprintf('Using %d instances of the reference RDM for random-effects test of RDM relatedness.\n',nRefRDMinstances);

    elseif isequal(subjectRFXacross,'candRDMinstances')
        %         fprintf('Found less than 12 instances of the reference RDM.'); %
        %         we are already saying this earlier
        fprintf('Using %d instances of all candidate RDMs for random-effects  test of RDM relatedness.\n',nsCandRDMinstances(1));
    else
        userOptions.RDMrelatednessTest='randomisation';
    end
end

switch userOptions.RDMrelatednessTest,
    case 'subjectRFXsignedRank'
        fprintf('\nPerforming signed-rank test of RDM relatedness (subject as random effect).\n');
        for candI = 1:numel(candRDMs)
            [ps(candI)] = signrank_onesided(cand2refSims(:,candI));
        end
        ps = ps(sortedIs);
        stats_p_r.candRelatedness_p = ps;
        nSubjects = size(cand2refSims,1);
        es = std(cand2refSims)/sqrt(nSubjects);
        es = es(sortedIs);



    case 'subjectRFXbootstrap'
        fprintf('\nPerforming subject bootstrap test of RDM relatedness (subject as random effect).\n');
        userOptions.resampleSubjects=1;userOptions.resampleConditions=0;

        [realRs bootstrapEs pairwisePs_bootSubj bootstrapRs] = ...
            bootstrapRDMs( ...
            refRDM_stack, meanCandRDMs, userOptions);

        % bootstrapRs is now nCandRDMs x nBootstrap
        bootstrapRs=bootstrapRs';% transpose to make it nBootstrap
                                 % x nCandRDMs
        for candI = 1:numel(candRDMs)
            [ci_lo(candI), ci_up(candI), ps(candI)] = bootConfint( ...
                cand2refSims(:,candI), bootstrapRs(:,candI), 'greater', ...
                userOptions);
        end
        ps = ps(sortedIs);
        stats_p_r.candRelatedness_p = ps;
        es = std(bootstrapRs);
        es = es(sortedIs);


    case 'subjectConditionRFXbootstrap'
        fprintf('\nPerforming subject and condition bootstrap test of RDM relatedness (subject and condition as random effects).\n');
        userOptions.resampleSubjects=1;userOptions.resampleConditions=1;

        [realRs bootstrapEs pairwisePs_bootSubj bootstrapRs] = ...
            bootstrapRDMs( ...
            refRDM_stack, meanCandRDMs, userOptions);

        % bootstrapRs is now nCandRDMs x nBootstrap
        bootstrapRs=bootstrapRs';% nBootstrap x nCandRDMs
        for candI = 1:numel(candRDMs)
            [ci_lo(candI), ci_up(candI), ps(candI)] = bootConfint( ...
                cand2refSims(:,candI), bootstrapRs(:,candI), 'greater', ...
                userOptions);
        end
        ps = ps(sortedIs);
        stats_p_r.candRelatedness_p = ps;
        es = std(bootstrapRs);es = es(sortedIs);


    case 'randomisation'
        fprintf('\nPerforming condition-label randomisation test of RDM relatedness (fixed effects).\n');
        nRandomisations = userOptions.nRandomisations;
        for rdmI = 1:numel(candRDMs)
            % do the randomisation test, also keep the randomistion correltions
            % in a separte matrix
            rdms(rdmI,:) = vectorizeRDM(meanCandRDMs(:,:,rdmI));
        end
        [n,n]=size(meanRefRDM);
        exhaustPermutations = false;
        if n < 8
            allPermutations = exhaustivePermutations(n);
            nRandomisations = size(allPermutations, 1);
            exhaustPermutations = true;
            warning('(!) Comparing RDMs with fewer than 8 conditions (per conditions set) will produce unrealiable results!\n  + I''ll partially compensate by using exhaustive instead of random permutations...');
        end%if
        % make space for null-distribution of correlations
        rs_null=nan(userOptions.nRandomisations,numel(candRDMs));

        % index method would require on the order of n^2*nRandomisations
        % memory, so i'll go slowly for now...
        %tic
        if isequal(userOptions.RDMcorrelationType,'Kendall_taua')
            for randomisationI=1:userOptions.nRandomisations
                if exhaustPermutations
                    randomIndexSeq = allPermutations(randomisationI, :);
                else
                    randomIndexSeq = randomPermutation(n);
                end%if

                rdmA_rand_vec=vectorizeRDM(meanRefRDM(randomIndexSeq,randomIndexSeq));
                for candI = 1:numel(candRDMs)
                    rs_null(randomisationI,candI)=rankCorr_Kendall_taua(rdmA_rand_vec',rdms(candI,:)');
                end
                if mod(randomisationI,floor(userOptions.nRandomisations/100))==0
                    fprintf('%d%% ',floor(100*randomisationI/userOptions.nRandomisations))
                    if mod(randomisationI,floor(userOptions.nRandomisations/10))==0, fprintf('\n'); end;
                end
            end % randomisationI
            fprintf('\n');
        elseif isequal(userOptions.RDMcorrelationType,'raeSpearman')
            for randomisationI=1:userOptions.nRandomisations
                if exhaustPermutations
                    randomIndexSeq = allPermutations(randomisationI, :);
                else
                    randomIndexSeq = randomPermutation(n);
                end%if

                rdmA_rand_vec=vectorizeRDM(meanRefRDM(randomIndexSeq,randomIndexSeq));
                for candI = 1:numel(candRDMs)
                    rs_null(randomisationI,candI)=raeSpearmanCorr(rdmA_rand_vec',rdms(candI,:)');
                end
                if mod(randomisationI,floor(userOptions.nRandomisations/100))==0
                    fprintf('%d%% ',floor(100*randomisationI/userOptions.nRandomisations))
                    if mod(randomisationI,floor(userOptions.nRandomisations/10))==0, fprintf('\n'); end;
                end
            end % randomisationI
            fprintf('\n');
        else
            for randomisationI=1:userOptions.nRandomisations
                if exhaustPermutations
                    randomIndexSeq = allPermutations(randomisationI, :);
                else
                    randomIndexSeq = randomPermutation(n);
                end%if

                rdmA_rand_vec=vectorizeRDM(meanRefRDM(randomIndexSeq,randomIndexSeq));
                rs_null(randomisationI,:)=corr(rdmA_rand_vec',rdms','type',userOptions.RDMcorrelationType,'rows','pairwise');
                if mod(randomisationI,floor(userOptions.nRandomisations/100))==0
                    fprintf('%d%% ',floor(100*randomisationI/userOptions.nRandomisations))
                    if mod(randomisationI,floor(userOptions.nRandomisations/10))==0, fprintf('\n'); end;
                end
            end % randomisationI
            fprintf('\n');
        end

        % p-values from the randomisation test
        for candI = 1:numel(candRDMs)
            p_randCondLabels(candI) = 1 - relRankIn_includeValue_lowerBound(rs_null(:,candI),y(candI)); % conservative
        end
        p_randCondLabels = p_randCondLabels(sortedIs);
        % p-values corrected for MC by controlling FWE rate
        for candI = 1:numel(candRDMs)
            p_randCondLabels_fwe(candI) = 1 - relRankIn_includeValue_lowerBound(max(rs_null'),y(candI)); % conservative
        end
        p_randCondLabels_fwe = p_randCondLabels_fwe(sortedIs);
        stats_p_r.candRelatedness_p = [p_randCondLabels;p_randCondLabels_fwe];

    otherwise
        fprintf('Not performing any test of RDM relatedness.\n');
end


%% decide inference type and perform candidate-RDM-comparison test
if isequal(userOptions.candRDMdifferencesTest,'subjectRFXsignedRank') || isequal(userOptions.candRDMdifferencesTest,'subjectRFXbootstrap') || isequal(userOptions.RDMrelatednessTest,'subjectConditionRFXbootstrap')
    % user requested subject random-effects inference (by signed-rank or bootstrap test)
    switch subjectRFXacross
        case 'refRDMinstances';
            fprintf('Using %d instances of the reference RDM for random-effects tests comparing pairs of candidate RDMs.\n',nRefRDMinstances);

        case 'candRDMinstances'
            fprintf('Using %d instances of all candidate RDMs for random-effects tests comparing pairs of candidate RDMs.\n',nsCandRDMinstances(1));

        otherwise
            fprintf('Attempting to revert to condition-bootstrap tests comparing pairs of candidate RDMs.\n');
            if conditionRFX
                fprintf('reverting to condition bootstrap tests for comparing pairs of candidate RDMs\n');
                userOptions.candRDMdifferencesTest='conditionRFXbootstrap';
            else
                fprintf('there are less than 20 conditions. can not do tests for comparing pairs of candidate RDMs\n');
                userOptions.candRDMdifferencesTest='none';
            end
    end
end

switch userOptions.candRDMdifferencesTest,
    case 'subjectRFXsignedRank'
        fprintf('Performing signed-rank test for candidate RDM comparisons (subject as random effect).\n');
        % do a one sided signrank test on the similarity of the each candidate
        % RDM,averaged over all instances with the different instances of the
        % reference RDM
        for candI1 = 1:numel(candRDMs)
            for candI2 = 1:numel(candRDMs)
                [pairWisePs(candI1,candI2)] = signrank(cand2refSims(:,candI1),cand2refSims(:,candI2),'alpha',0.05,'method','exact');
            end
        end
        stats_p_r.candDifferences_p = pairWisePs(sortedIs,sortedIs);
        stats_p_r.SEs=es;

        hold on
        errorbar(1:numel(candRDMs),y_sorted,es,'Color',[0 0 0],'LineWidth',...
                2,'LineStyle','none');

    case 'subjectRFXbootstrap'
        fprintf('Performing subject bootstrap test comparing candidate RDMs (subject as random effect).\n');
        if ~exist('pairwisePs_bootSubj','var')
            [realRs bootstrapEs pairwisePs_bootSubj bootstrapRs] = bootstrapRDMs(refRDM_stack, meanCandRDMs, userOptions);
            for candI = 1:numel(candRDMs)
                [ci_lo(candI), ci_up(candI)] = bootConfint( ...
                    cand2refSims(:,candI), bootstrapRs(:,candI), 'greater', ...
                    userOptions);
            end
        end
        hold on
        errorbar(1:numel(candRDMs),y_sorted, y_sorted - ci_lo(sortedIs), ...
            ci_up(sortedIs) - y_sorted,'Color',[0 0 0],'LineWidth',...
            2,'LineStyle','none');

        pairWisePs = nan(numel(candRDMs));
        for candI = 1:(numel(candRDMs)-1)
            bsI = bootstrapRs(:,candI);
            cand2refSimsI = cand2refSims(:,candI);
            for candJ = (candI+1):numel(candRDMs)
                bsJ = bootstrapRs(:,candJ);
                cand2refSimsJ = cand2refSims(:,candJ);
                [~, ~, pairWisePs(candI, candJ)] = bootConfint( ...
                    cand2refSimsI - cand2refSimsJ, bsI - bsJ, ...
                    'two-tailed', userOptions);
                pairWisePs(candJ, candI) = pairWisePs(candI, candJ);
            end
        end
        stats_p_r.candDifferences_p = pairWisePs(sortedIs,sortedIs);
        stats_p_r.SEs = bootstrapEs(sortedIs);

    case 'subjectConditionRFXbootstrap'
        fprintf('Performing condition and subject bootstrap test comparing candidate RDMs (condition and subject as random effects).\n');
        if ~exist('pairwisePs_bootSubjCond','var')
            [realRs bootstrapEs pairwisePs_bootSubjCond bootstrapRs] = bootstrapRDMs(refRDM_stack, meanCandRDMs, userOptions);
        end

        hold on
        errorbar(1:numel(candRDMs),y_sorted,bootstrapEs(sortedIs),'Color',[0 0 0],'LineWidth',...
                2,'LineStyle','none');

        pairWisePs = pairwisePs_bootSubjCond;
        stats_p_r.candDifferences_p = pairWisePs(sortedIs,sortedIs);
        stats_p_r.SEs = bootstrapEs(sortedIs);

    case 'conditionRFXbootstrap'
        fprintf('Performing condition bootstrap test comparing candidate RDMs (subject as random effect).\n');
        userOptions.resampleSubjects=0;userOptions.resampleConditions=1;
        [realRs bootstrapEs pairwisePs bootstrapRs] = bootstrapRDMs(refRDM_stack, meanCandRDMs, userOptions);
        for candI = 1:numel(candRDMs)
            [ci_lo(candI), ci_up(candI)] = bootConfint( ...
                cand2refSims(:,candI), bootstrapRs(:,candI), 'greater', ...
                userOptions);
        end
        pairWisePs = nan(numel(candRDMs));
        for candI = 1:(numel(candRDMs)-1)
            bsI = bootstrapRs(:,candI);
            cand2refSimsI = cand2refSims(:,candI);
            for candJ = (candI+1):numel(candRDMs)
                bsJ = bootstrapRs(:,candJ);
                cand2refSimsJ = cand2refSims(:,candJ);
                [~, ~, pairWisePs(candI, candJ)] = bootConfint( ...
                    cand2refSimsI - cand2refSimsJ, bsI - bsJ, ...
                    'two-tailed', userOptions);
                pairWisePs(candJ, candI) = pairWisePs(candI, candJ);
            end
        end
        stats_p_r.candDifferences_p = pairwisePs(sortedIs,sortedIs);
        stats_p_r.SEs = bootstrapEs(sortedIs);
        hold on
        errorbar(1:numel(candRDMs),y_sorted,bootstrapEs(sortedIs),'Color',[0 0 0],'LineWidth',...
                2,'LineStyle','none');

    otherwise
        fprintf('Not performing any test for comparing pairs of candidate RDMs.\n');
end

if min(y_sorted)<0
    if ~exist('es','var')
        es = 0.1*ones(size(y_sorted));
    end
    Ymin=min(y_sorted)-max(es)-.05;
else
    Ymin=0;
end


%% display the ceiling
set(gcf,'Renderer','OpenGL');

h=patch([0.1 0.1 numel(candRDMs)+1 numel(candRDMs)+1],[ceilingLowerBound ceilingUpperBound ceilingUpperBound ceilingLowerBound],[0.1 0.1 0.1 0.1],[.7 .7 .7],'edgecolor','none');hold on;
alpha(h,0.5);


%% add the p-values from RDM relatedness tests
if ~isequal(userOptions.RDMrelatednessTest,'none')
    ps = stats_p_r.candRelatedness_p(1,:);

    pStringCell = cell(1, numel(candRDMs));

    if isequal(userOptions.RDMrelatednessMultipleTesting,'none') || isequal(userOptions.RDMrelatednessMultipleTesting,'FWE')
        thresh = userOptions.RDMrelatednessThreshold;
    elseif isequal(userOptions.RDMrelatednessMultipleTesting,'FDR')
        thresh = FDRthreshold(ps,userOptions.candRDMdifferencesThreshold);
    end

    if isequal(userOptions.RDMrelatednessMultipleTesting,'FWE')
        ps_fwe = stats_p_r.candRelatedness_p(2,sortedIs);

        for test = 1:numel(candRDMs)
            if isequal(userOptions.plotpValues,'*')
                if ps(test) <= ps_fwe(test)
                    pStringCell{test} = '*';
                else
                    pStringCell{test} = '';
                end
            else
                if test==1
                    pStringCell{test} = ['p = ' num2str(ps(test), 2)]; % 2 significant figures!
                else
                    pStringCell{test} = [num2str(ps(test), 2)]; % 2 significant figures!
                end
            end
        end
    else
        for test = 1:numel(candRDMs)
            if isequal(userOptions.plotpValues,'*')
                if ps(test) <= thresh
                    pStringCell{test} = '*';
                else
                    pStringCell{test} = '';
                end
            else
                if test==1
                    if ps(test) <= thresh
                        pStringCell{test} = ['\bf' compactPvalueString(ps(test))];
                    else
                        pStringCell{test} = ['\rmp = ' compactPvalueString(ps(test))];
                    end
                else
                    pStringCell{test} = [compactPvalueString(ps(test))]; % 2 significant figures!
                end
                if ps(test) <= thresh
                    pStringCell{test} = ['\bf',pStringCell{test}];
                else
                    pStringCell{test} = ['\rm',pStringCell{test}];
                end
            end%switch:pIndication
        end%for:test
    end % isequal(userOptions.RDMrelatednessTest....

    for test = 1:numel(candRDMs)
        switch userOptions.plotpValues
            case '*'
                text(test, -0.06, ['\bf\fontsize{14}',deunderscore(pStringCell{test})], 'Rotation', 0, 'Color', [0 0 0],'HorizontalAlignment','center');
            case '='
                text(test, -0.06, ['\fontsize{8}',deunderscore(pStringCell{test})], 'Rotation', 0, 'Color', [0 0 0],'HorizontalAlignment','center');
        end
    end
end
if isequal(userOptions.plotpValues,'=')
    text(.5,-0.06,['\fontsize{8}\rmp= '],'Rotation', 0, 'Color', [0 0 0],'HorizontalAlignment','center')
end


%% label the bars with the names of the candidate RDMs
for test = 1:numel(candRDMs)
    text(test, Ymin-0.05, ['\bf',deunderscore(candidateNames{test})], 'Rotation', 45, 'Color', [0 0 0],'HorizontalAlignment','right');
end


%% add the pairwise comparison horizontal lines
thresh = userOptions.candRDMdifferencesThreshold;
switch userOptions.candRDMdifferencesMultipleTesting
    case 'FDR'
        pMat = stats_p_r.candDifferences_p;
        pMat(logical(eye(numel(candRDMs)))) = 0;
        allPairwisePs = squareform(pMat);
        threshold = FDRthreshold(allPairwisePs,thresh);

    case 'FWE'
        nTests = numel(candRDMs)*(numel(candRDMs)-1)/2;
        threshold = thresh/nTests;
    case 'none'
        threshold = thresh;
end
if nSubjects == 1
    yy=addComparisonBars(stats_p_r.candDifferences_p,(max(y)+0.1),threshold);
else
    yy=addComparisonBars(stats_p_r.candDifferences_p,(ceilingUpperBound+0.1),threshold);
end
cYLim = get(gca, 'YLim');
cYMax = cYLim(2);
labelBase = cYMax + 0.1;
nYMax = labelBase;
nYLim = [0, nYMax];
set(gca, 'YLim', nYLim);
hold on;

maxYTickI=ceil(max([y_sorted ceilingUpperBound])*10);
% set(gca,'YTick',[0 1:maxYTickI]./10);
% set(gca,'XTick',[]);
axis off;

% plot pretty vertical axis
lw=1;
for YTickI=0:maxYTickI
    plot([0 0.1],[YTickI YTickI]./10,'k','LineWidth',lw);
    text(0,double(YTickI/10),num2str(YTickI/10,1),'HorizontalAlignment','right');
end
plot([0.1 0.1],[0 YTickI]./10,'k','LineWidth',lw);
text(-1.2,double(maxYTickI/10/2),{'\bf RDM correlation ',['\rm[',deunderscore(userOptions.RDMcorrelationType),', averaged across ',num2str(nRefRDMinstances),' subjects]']},'HorizontalAlignment','center','Rotation',90);

if ~exist('yy','var')
    ylim([Ymin max(y_sorted)+max(es)+.05]);
else
    ylim([Ymin yy+.1]);
end

if exist('refRDMName','var') & nRefRDMinstances == 1
    title({['\bf\fontsize{11} How closely is ',refRDMName,' related to each of the candidate RDMs?'],['\fontsize{9}RDM relatedness tests: \rm',userOptions.RDMrelatednessTest,...
    '\rm (threshold: ',num2str(userOptions.RDMrelatednessThreshold),', multiple testing: ',userOptions.RDMrelatednessMultipleTesting,')'],...
    ['\bfpairwise comparison tests: \rm',userOptions.candRDMdifferencesTest,...
    '\rm (threshold: ',num2str(userOptions.candRDMdifferencesThreshold),', multiple testing: ',userOptions.candRDMdifferencesMultipleTesting,')']})
else
    title({'\bf\fontsize{11} How closely is the reference RDM related to each of the candidate RDMs?',['\fontsize{9}RDM relatedness tests: \rm',userOptions.RDMrelatednessTest,...
    '\rm (threshold: ',num2str(userOptions.RDMrelatednessThreshold),', multiple testing: ',userOptions.RDMrelatednessMultipleTesting,')'],...
    ['\bfpairwise comparison tests: \rm',userOptions.candRDMdifferencesTest,...
    '\rm (threshold: ',num2str(userOptions.candRDMdifferencesThreshold),', multiple testing: ',userOptions.candRDMdifferencesMultipleTesting,')']})
end
axis square
h=get(gcf,'Position');
h(4)=h(4)+100;

axis([0 numel(candRDMs)+0.5 -0.5 yy+0.1]);
%rectangle('Position',[0 -0.5 numel(candRDMs)+0.5 0.5+yy+0.1]);
pageFigure(userOptions.figureIndex(1));
handleCurrentFigure([userOptions.resultsPath,filesep,userOptions.figure1filename],userOptions);


%% display the matrices for pairwise candidate RDM comparisons
pMask = true(size(stats_p_r.candDifferences_r));
pMask(stats_p_r.candDifferences_r > 0) = false;% the entries that are '1' will not be displayed
h = figure(userOptions.figureIndex(2)); clf;
set(h, 'Color', 'w');
clear pMat; pMat = stats_p_r.candDifferences_p;
% display the uncorrected p-values in the first figure
nCols = 256;clims = [0 1];
cols = colorScale([1 0 0;0 0 0],nCols);
pMat(pMask)= 1;% arbitrary large number, would be displayed as black
subplot(221)
RGBim=mat2RGBimage(pMat,cols,clims);
for k1 = 1:size(pMat,1)
    for k2 = 1:size(pMat,2)
        if pMask(k1,k2)
            RGBim(k1,k2,:) = [1 1 1];
        end
    end
end
image(RGBim);
colormap(cols);axis square
set(gca,'xTick',1:numel(candRDMs),'xTickLabel',candidateNames(1:numel(candRDMs)));
xticklabel_rotate
set(gca,'yTick',1:numel(candRDMs),'yTickLabel',candidateNames(1:numel(candRDMs)));
title({'\bfuncorrected p-values','\rmred-to-black colors linearly map to the p-value range [0 1]'});

% uncorrected thresholded pMat
thr = userOptions.candRDMdifferencesThreshold;
u = pMat;
u(pMat>thr) = 1;% black
u(pMat <= thr/5) = 0;% red
u(pMat<=thr & pMat>thr/5) = 0.5;% intermediate color
subplot(222)
RGBim=mat2RGBimage(u,cols,clims);
for k1 = 1:size(pMat,1)
    for k2 = 1:size(pMat,2)
        if pMask(k1,k2)
            RGBim(k1,k2,:) = [1 1 1];
        end
    end
end
image(RGBim);
colormap(cols);axis square
set(gca,'xTick',1:numel(candRDMs),'xTickLabel',candidateNames(1:numel(candRDMs)));
xticklabel_rotate
set(gca,'yTick',1:numel(candRDMs),'yTickLabel',candidateNames(1:numel(candRDMs)));
title({'\bfthresholded uncorrected p-values',['\rmdark: n.s , dark red: p<',num2str(thr,'%0.3f'),', red: p<',num2str(thr/5)]});

% FDR control
Ps = stats_p_r.candDifferences_p;
Ps(logical(eye(numel(candRDMs)))) = 0;
allPairwisePs = squareform(Ps);
thr_fdr = FDRthreshold(allPairwisePs,userOptions.candRDMdifferencesThreshold);
thr_fdr2 = FDRthreshold(allPairwisePs,userOptions.candRDMdifferencesThreshold/5);
u = pMat;
u(pMat>thr_fdr) = 1;% black
u(pMat <= thr_fdr/5) = 0;% red
u(pMat<=thr_fdr & pMat>thr_fdr/5) = 0.5;% intermediate color
subplot(223)
RGBim=mat2RGBimage(u,cols,clims);
for k1 = 1:size(pMat,1)
    for k2 = 1:size(pMat,2)
        if pMask(k1,k2)
            RGBim(k1,k2,:) = [1 1 1];
        end
    end
end
image(RGBim);
colormap(cols);axis square
set(gca,'xTick',1:numel(candRDMs),'xTickLabel',candidateNames(1:numel(candRDMs)));
xticklabel_rotate
set(gca,'yTick',1:numel(candRDMs),'yTickLabel',candidateNames(1:numel(candRDMs)));
title({'\bfthresholded to control the FDR',['\rmdark: n.s., dark red: E(FDR)<',num2str(userOptions.candRDMdifferencesThreshold,'%0.3f'),', red: E(FDR)<',num2str(userOptions.candRDMdifferencesThreshold/5,'%0.4f')]});

% FWE (Bonferroni) control
thr_bnf = 2*userOptions.candRDMdifferencesThreshold/(numel(candRDMs)*(numel(candRDMs)-1));
u = pMat;
u(pMat>thr_bnf) = 1;% black
u(pMat <= thr_bnf/5) = 0;% red
u(pMat<=thr_bnf & pMat>thr_bnf/5) = 0.5;% intermediate color
subplot(224)
RGBim=mat2RGBimage(u,cols,clims);
for k1 = 1:size(pMat,1)
    for k2 = 1:size(pMat,2)
        if pMask(k1,k2)
            RGBim(k1,k2,:) = [1 1 1];
        end
    end
end
image(RGBim);
colormap(cols);axis square
set(gca,'xTick',1:numel(candRDMs),'xTickLabel',candidateNames(1:numel(candRDMs)));
xticklabel_rotate
set(gca,'yTick',1:numel(candRDMs),'yTickLabel',candidateNames(1:numel(candRDMs)));
title({'\bfthresholded to control the FWE (Bonferroni)',['\rmblack: n.s., dark red: p(corr.)<',num2str(userOptions.candRDMdifferencesThreshold,'%0.4f'),', red: p(corr.)<',num2str(userOptions.candRDMdifferencesThreshold/5,'%0.4f')]});

handleCurrentFigure([userOptions.resultsPath,filesep,userOptions.figure2filename],userOptions);


%% save it
cd(fullfile(userOptions.rootPath));
if userOptions.saveFiguresPS
    saveas(gcf,[userOptions.analysisName,'_comparingRefRDM2CandRDMs','.eps'],'eps');
end
if userOptions.saveFiguresFig
    saveas(gcf,[userOptions.analysisName,'_comparingRefRDM2CandRDMs','.fig'],'fig');
end
if userOptions.saveFiguresPDF
    exportCurrentFigAsPDF([userOptions.analysisName,'_comparingRefRDM2CandRDMs'],userOptions);
end
if userOptions.saveFiguresPS
    exportCurrentFigAsPostscript([userOptions.analysisName,'.ps'],userOptions);
end

end%function
