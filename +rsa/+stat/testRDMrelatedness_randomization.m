function [p_randCondLabels,r,p_randDissims,p_conv]=testRDMrelatedness_randomization(rdmA,rdmB,conditionSetIs_vector,corrType)

% USAGE
%       [p_randCondLabels,r,p_randDissims,p_conv]=
%       testRDMRelatedness_rankCorr_randomization(rdmA,rdmB
%       [,conditionSetsVector_LOG])
%
% FUNCTION
%       tests the null hypothesis that similarity matrices A and B are
%       unrelated. the test simulates a null distribution of correlations
%       between A and B by means of randomization of the conditions labels
%       of the similarity matrix.


%% control variables
nRandomizations=10000; % 10,000


%% preparations
% square the RDMs
rdmA=stripNsquareRDMs(rdmA);
rdmB=stripNsquareRDMs(rdmB);
[n,n]=size(rdmA);
[nB,nB]=size(rdmB);
if n~=nB
    error('testRDMRelatedness_randomization: RDMs need to be of the same size.');
end

if ~exist('conditionSetIs_vector','var'), conditionSetIs_vector=ones(1,n); end
conditionSetIndices=unique(conditionSetIs_vector(conditionSetIs_vector~=0)); % '0' indicates: to be excluded
nConditionSets=numel(conditionSetIndices);

if ~exist('corrType','var'), corrType='Spearman'; end;


%% test relatedness within a single conditions set
if nConditionSets==1;

    % reduce the RDMs
    conditionSet_LOG=(conditionSetIs_vector==conditionSetIndices(1));
    rdmA=rdmA(conditionSet_LOG,conditionSet_LOG);
    rdmB=rdmB(conditionSet_LOG,conditionSet_LOG);
    [n,n]=size(rdmA);
    
    % vectorize the RDMs
    rdmA_vec=vectorizeRDM(rdmA);
    rdmB_vec=vectorizeRDM(rdmB);
  
    % make space for null-distribution of correlations
    rs=nan(nRandomizations,1);
    
    % index method would require on the order of n^2*nRandomizations
    % memory, so i'll go slowly for now...
    %tic
    for randomizationI=1:nRandomizations
        randomIndexSeq=randperm(n);
        rdmA_rand_vec=vectorizeRDM(rdmA(randomIndexSeq,randomIndexSeq));
        
        rs(randomizationI)=corr(rdmA_rand_vec',rdmB_vec','type',corrType,'rows','pairwise');  % correlation types: Pearson (linear), Spearman (rank), Kendall (rank)
        % setting the 'rows' parameter to 'pairwise' makes the function ignore nans.
    end
    %toc

end


%% test relatedness between two disjoint conditions sets
if nConditionSets==2;

    % reduce the RDMs (they become nonsquare: different input and output sets)
    conditionSet1_LOG=(conditionSetIs_vector==conditionSetIndices(1));
    conditionSet2_LOG=(conditionSetIs_vector==conditionSetIndices(2));

    rdmA=rdmA(conditionSet1_LOG,conditionSet2_LOG);
    rdmB=rdmB(conditionSet1_LOG,conditionSet2_LOG);
    [n1,n2]=size(rdmA);
    
    % vectorize the RDMs (already no diagonal-including sections present)
    rdmA_vec=rdmA(:)';
    rdmB_vec=rdmB(:)';
  
    % make space for null-distribution of correlations
    rs=nan(nRandomizations,1);
    
    % index method would require on the order of n^2*nRandomizations
    % memory, so i'll go slowly for now...
    %tic
    for randomizationI=1:nRandomizations
        conditionsSet1_randomIndexSeq=randperm(n1);
        conditionsSet2_randomIndexSeq=randperm(n2);
        rdmA_rand=rdmA(conditionsSet1_randomIndexSeq,conditionsSet2_randomIndexSeq);
        rdmA_rand_vec=rdmA_rand(:)';
        
        rs(randomizationI)=corr(rdmA_rand_vec',rdmB_vec','type',corrType,'rows','pairwise');  % correlation types: Pearson (linear), Spearman (rank), Kendall (rank)
        % setting the 'rows' parameter to 'pairwise' makes the function ignore nans.
    end
    %toc

end


%% compute actual correlation and the corresponding conventional p value
[r, p_conv]=corr(rdmA_vec',rdmB_vec','type',corrType,'rows','pairwise');  % correlation types: Pearson (linear), Spearman (rank), Kendall (rank)
% setting the 'rows' parameter to 'pairwise' makes the function ignore nans.


%% compute p value for condition-label randomization
% this is a valid method.
p_randCondLabels=1-relRankIn_includeValue_lowerBound(rs,r); % conservative


%% compute p value for dissimilarity-set randomization
% valid method???
% nDissimilarities=numel(rdmA_vec);
% rs=nan(nRandomizations,1);
% 
% for randomizationI=1:nRandomizations
%     rdmA_rand_vec=rdmA_vec(randperm(nDissimilarities));
% 
%     rs(randomizationI)=corr(rdmA_rand_vec',rdmB_vec','type','Spearman','rows','pairwise');  % correlation types: Pearson (linear), Spearman (rank), Kendall (rank)
%     % setting the 'rows' parameter to 'pairwise' makes the function ignore nans.
% end
% 
% p_randDissims=1-relRankIn_includeValue_lowerBound(rs,r); % conservative

p_randDissims=nan;

