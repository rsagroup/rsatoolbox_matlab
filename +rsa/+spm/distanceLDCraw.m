function [d, Sig] = distanceLDCraw(Y,SPM,conditionVec,varargin)
% function d=rsa.spm.distanceLDCraw(Y,SPM,conditionVec,varargin);
% First, gets the regression coefficent from the SPM, and prewhitens them.
% Prewhiten can be controlled using different methods (run-wise or overall).
% The same code is used in rsa.spm.noiseNormalizeBeta.
% Secondly, it calcualtes LDC stances with a leave-one out crossvalidation,
% It uses optimal combination of the beta-coeeficients in the part that
% averages across partitions. By default, the different partitions are
% assumed to be the different imaging runs.
% Note that in calculating the distances, the structure of the first-level 
% design matrix is optimally taken into account. The command: 
%  RDM  = rsa.spm.distanceLDCraw(Y,SPM,condition); 
% is therefore equivalent to:  
%  beta = rsa.spm.noiseNormalizeBeta(Y,SPM); 
%  RDM  = rsa.distanceLDC(beta,partition,condition,SPM.xX.xKXs.X); 
% INPUT:
%    Y             raw timeseries, T by P
%    SPM:          SPM structure
%    conditionVec: Vector that indicates conditions for each column of design matrix 
%                  set to 0 for nuisance regressors
% OUTPUT:
%    d:            numCond*(numCond-1)/2 distances between experimental
%                  conditions
% OPTIONs:
%   'normmode':    'runwise': Does the multivariate noise normalisation by run
%                  'overall': Does the multivariate noise normalisation overall
%   'normmethod':  'multivariate': The is the default using ledoit-wolf reg.
%                  'univariate': Performing univariate noise normalisation (t-values)
%                  'none':    No noise normalisation
% (c) 2015 Joern Diedrichsen, Alex Walther
Opt.normmode = 'runwise';  % Either runwise or overall
Opt.normmethod = 'multivariate';  % Either runwise or overall
Opt = rsa.getUserOptions(varargin,Opt,{'normmode','normmethod'});

[T,numVox]=size(Y);                                             %%% number of time points and voxels

%%% Discard NaN voxels
test=isnan(sum(Y));
if (any(test))
    warning(sprintf('%d of %d voxels contained NaNs -discarding',sum(test),length(test)));
    Y=Y(:,test==0);
end;

xX    = SPM.xX;                                            %%% take the design
X     = SPM.xX.xKXs.X;
numReg = size(X,2);

% Check condition vector
numCond = max(conditionVec);
if (length(conditionVec)<numReg)
    conditionVec=[conditionVec;zeros(numReg-length(conditionVec),1)];
end;
Z = indicatorMatrix('identity_p',conditionVec);
nonInterest = all(Z==0,2);   % Regressors not in the conditions
numNonInterest = sum(nonInterest);
Z(nonInterest,end+1:end+sum(numNonInterest))=eye(numNonInterest);
C = indicatorMatrix('allpairs',[1:numCond]);

%%% Get partions: For each run (1:K), find the time points (T) and regressors (K+Q) that belong to the run
partT = nan(T,1);
partN = nan(numReg,1);
numPart=length(SPM.Sess);                                     %%% number of runs
for i=1:numPart
    partT(SPM.Sess(i).row,1)=i;
    partN(SPM.Sess(i).col,1)=i;
    partN(SPM.xX.iB(i),1)=i;                                %%% Add intercepts
end;

KWY=spm_filter(xX.K,xX.W*Y);                               %%% filter out low-frequence trends in Y
res=spm_sp('r',xX.xKXs,KWY);                               %%% residuals: res  = Y - X*beta

%%% do run-wise noise normalization
switch (Opt.normmethod)
    case 'none'
        % Do nothing 
    case 'multivariate'
        switch (Opt.normmode)
            case 'runwise'
                for i=1:numPart
                    idxT = partT==i;
                    idxN = partN==i;
                    numFilt = size(xX.K(i).X0,2);
                    sum(idxT)-sum(idxN)-numFilt-1;
                    [Sw_hat(:,:,i),shrink(i)]=covdiag(res(idxT,:),sum(idxT)-sum(idxN)-numFilt-1);                    %%% regularize Sw_hat through optimal shrinkage
                    [V,L]=eig(Sw_hat(:,:,i));       % This is overall faster and numerical more stable than Sw_hat.^-1/2
                    l=diag(L);
                    sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
                    KWY(idxT,:)=KWY(idxT,:)*sq;
                end;
            case 'overall'
                Sw_hat=covdiag(res,SPM.xX.erdf);    % regularize Sw_hat through optimal shrinkage
                [V,L]=eig(Sw_hat);                  % This is overall faster and numerical more stable than Sw_hat.^-1/2
                l=diag(L);
                sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
                KWY=KWY*sq;
        end;
    case 'univariate'
        error('univariate noise normalisation not implemented yet');
    otherwise
        error('normmethod needs to be ''multivariate'', ''univariate'', or ''none''');
end;

% Estimate condition means within each
A = zeros(length(nonInterest(partN==partN(1))),numVox,numPart);
for i=1:numPart
    % Get the betas from the test run 
    indxN = partN==i;
    indxT = partT==i;
    Za = Z(indxN,:);
    Za = Za(:,any(Za,1));
    Xa = X(indxT,indxN);
    Ma  = Xa*Za;
    A(:,:,i)     = (Ma'*Ma)\(Ma'*KWY(indxT,:));
    
    % Get the betas based on the other runs 
    indxN = partN~=i;
    indxT = partT~=i;
    Zb    = Z(indxN,:);
    Zb    = Zb(:,any(Zb,1));
    Xb    = X(indxT,indxN);
    Mb    = Xb*Zb;
    B     = (Mb'*Mb)\Mb'*KWY(indxT,:);
    
    % Pick condition of interest
    interest = ~nonInterest(partN==i);
    
    % Caluclate distances 
    %d(i,:)= sum((C*A(1:numCond,:,i)).*(C*B(1:numCond,:)),2)'/numVox;      % Note that this is normalised to the number of voxels
    d(i,:)= sum((C*A(interest,:,i)).*(C*B(interest,:)),2)'/numVox;      % Note that this is normalised to the number of voxels
end;
d = sum(d)./numPart;

% If requested, also calculate the estimated variance-covariance 
% matrix from the residual across folds. 
A = A(interest,:,:);
if (nargout>1) 
    R=bsxfun(@minus,A,sum(A,3)/numPart);
    for i=1:numPart
        Sig(:,:,i)=R(:,:,i)*R(:,:,i)'/numVox;
    end;
    Sig=sum(Sig,3)/(numPart-1);
end; 
