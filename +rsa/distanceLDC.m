function [d,Sig]=distanceLDC(B,partition,conditionVec,X)
% function [d,Sigma]=distanceLDC(Y,partition,conditionVec,X);
% Calculates the crossvalidated squared Euclidian distances between activation patterns 
% (estimated beta weights) from the first level analysis. The distance is calculated 
% by using data across partitions (imaging runs), making these distance estimates unbiases 
% - i.e. the expected value is 0 if the two patterns are not different. This also 
% means that the distance estimates can be negative (see Walther et al., 2015).  
% If the betas weights are prewhitened, then the resultant distance will be a
% crossvalidated Mahalanobis distance (LDC). 
% If the optional input argument X is given, then the data will be combined
% across partitions in an optimal way, taking into account the different
% variabilities of the estimators. In this case, it is also allowed that
% not every regressor is used in every partition (although this is not recommended).

% The two lines:  
%  beta = rsa.spm.noiseNormalizeBeta(Y,SPM); 
%  RDM = rsa.distanceLDC(beta,partition,condition,SPM.xX.xKXs.X); 
% are equivalent to: 
%  RDM = rsa.spm.distanceLDCraw(Y,SPM,condition); 
% 
% INPUT:
%  B           : Noise-normalized activation patterns, a N x P matrix
%                If X is provided for optimal weighting of the regressors,
%                it is important also to submit all regressors of
%                no-interest (i.e. intercepts, etc), such that their value
%                can be taken into account 
%  partition   : N x 1 integer value that indicates the partition for crossvalidation (typically run number)
%                These need to be between 1...M, zeros are being ignored 
%                if the Partition vector is shorter than N, it is assumed that the
%                last regressors are intercept for the different runs-as is
%                usual in SPM.
%  conditionVec: N x 1 vector of conditions, zeros will be ignored as
%                regressors of no interest. If conditionVec is shorter than 
%                N, it is assumed that all remaining numbers are 0. 
%  X           : T x N Design matrix that is used to estimate from the first
%                level 
% OUTPUT: 
%   d          : Average distance values for each of the K*(K-1)/2 pairwise 
%               differences.  
%               (1xC), this measure is normalised to the numer of voxels 
%   Sig        : a KxK covariance matrix of the beta estimates across
%               different imaging runs. 
% Alexander Walther, Daan Wesselink, Joern Diedrichsen 
% 2/2016 

import rsa.util.*; 

[N,numVox]   = size(B); 
part    = unique(partition)';
part(part==0) = []; % Ignore the zero partitions 

numPart = numel(part);
numCond   = max(conditionVec); 

% Check on design matrix 
if (nargin>3 && ~isempty(X))  
    numReg     = size(X,2);             % Number of regressors in the first-level design matrix 
    if (numReg ~=N) 
        error('For optimal integration of beta weights, all N regressors (including no-interest) need to be submitted in Y'); 
    end; 
end; 

% Check length of partition vector  
missing = N-length(partition); 
if missing > 0 
    partition  = [partition;[1:missing]']; % Asssume that these are run intercepts 
end; 

% Check length of condition vector 
if (length(conditionVec)<N)
    conditionVec=[conditionVec;zeros(N-length(conditionVec),1)];
end; 

A = zeros(numCond,numVox,numPart);           % Allocate memory 
C = indicatorMatrix('allpairs',[1:numCond]); % Make contrats matrix 

% Make second-level design matrix, pulling through the regressors of no-interest 
Z = indicatorMatrix('identity_p',conditionVec); 
numNonInterest = sum(conditionVec==0);      % Number of no-interest regressors 
Z(conditionVec==0,end+[1:numNonInterest])=eye(numNonInterest); 

% Estimate condition means within each run and crossvalidate 
for i=1:numPart 
    % Left-out partition 
    indxA = partition==part(i);
    Za    = Z(indxA,:); 
    Za    = Za(:,any(Za,1));       % restrict to regressors that are not all 0
    Ba    = B(indxA,:);            % Get regression coefficients 

    % remainder of conditions 
    indxB = partition~=part(i);
    Zb    = Z(indxB,:); 
    Zb    = Zb(:,any(Zb,1));    % Restrict to regressors that are not all 0 
    Bb    = B(indxB,:);
    
    % Use design matrix if present to get GLS estimate 
   if (nargin>3 & ~isempty(X))
        Xa      = X(:,indxA);
        Xb      = X(:,indxB);
        indxX   = any(Xa,1);    % Restrict to regressors that are used in this partition
        Za      = Xa*Za; 
        Za      = Za(:,indxX); 
        Zb      = Xb*Zb; 
        Ba      = Xa(:,indxX)*Ba(indxX,:);
        Bb      = Xb*Bb; 
   end; 
    a     = pinv(Za)*Ba;
    b     = pinv(Zb)*Bb;
    A(:,:,i) = a(1:numCond,:); 
    d(i,:)= sum((C*A(:,:,i)).*(C*b(1:numCond,:)),2)'/numVox;      % Note that this is normalised to the number of voxels 
end; 
d = sum(d)./numPart; 

% If requested, also calculate the estimated variance-covariance 
% matrix from the residual across folds. 
if (nargout>1) 
    R=bsxfun(@minus,A,sum(A,3)/numPart);
    for i=1:numPart
        Sig(:,:,i)=R(:,:,i)*R(:,:,i)'/numVox;
    end;
    Sig=sum(Sig,3)/(numPart-1);
end; 

