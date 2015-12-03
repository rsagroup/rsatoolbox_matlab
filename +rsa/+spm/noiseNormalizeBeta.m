function [u_hat,resMS,Sw_raw,beta_hat,shrink]=noiseNormalizeBeta(Y,SPM,varargin)
% function [u_hat,Sw_hat,resMS,beta_hat]=rsa_noiseNormalizeBeta(Y,SPM,varargin)
% Estimates beta coefficiencts beta_hat and residuals from raw time series Y
% Estimates the true activity patterns u_hat by applying noise normalization to beta_hat
% INPUT:
%    Y        raw timeseries, T by P
%    SPM:     SPM structure
% OPTIONS:
%   'normmode':   'runwise': Does the multivariate noise normalisation by run
%                 'overall': Does the multivariate noise normalisation overall
%   'normmethod': 'multivariate': The is the default using ledoit-wolf reg.
%                 'univariate': Performing univariate noise normalisation (t-values)
%                 'none':    No noise normalisation
% OUTPUT:
%    u_hat    estimated true activity patterns (beta_hat after multivariate noise normalization),
%    resMS    residual mean-square - diagonal of the Var-cov matrix, 1 by P
%    Sw_raw   overall voxel error variance-covariance matrix (PxP), before regularisation 
%    beta_hat estimated raw regression coefficients, K*R by P
%    shrink   applied shrinkage factor 
% Alexander Walther, Joern Diedrichsen
% joern.diedrichsen@googlemail.com
% 2/2015
Opt.normmode = 'runwise';  % Either runwise or overall
Opt = rsa.getUserOptions(varargin,Opt);
[T,P]=size(Y);                                             %%% number of time points and voxels

%%% Discard NaN voxels
test=isnan(sum(Y));
if (any(test))
    warning(sprintf('%d of %d voxels contained NaNs -discarding',sum(test),length(test)));
    Y=Y(:,test==0);
end;

xX    = SPM.xX;                                            %%% take the design
[T,Q] = size(xX.X);

%%% Get partions: For each run (1:K), find the time points (T) and regressors (K+Q) that belong to the run
partT = nan(T,1);
partQ = nan(Q,1);
Nrun=length(SPM.Sess);                                     %%% number of runs
for i=1:Nrun
    partT(SPM.Sess(i).row,1)=i;
    partQ(SPM.Sess(i).col,1)=i;
    partQ(SPM.xX.iB(i),1)=i;                                %%% Add intercepts
end;

KWY=spm_filter(xX.K,xX.W*Y);                               %%% filter out low-frequence trends in Y
beta_hat=xX.pKX*KWY;                                       %%% ordinary least squares estimate of beta_hat = inv(X'*X)*X'*Y
res=spm_sp('r',xX.xKXs,KWY);                               %%% residuals: res  = Y - X*beta
clear KWY XZ                                               %%% clear to save memory

switch (Opt.normmode)
    case 'runwise'
        u_hat   = zeros(size(beta_hat));
        % do run-wise noise normalization
        shrink=zeros(Nrun,1);
        for i=1:Nrun
            idxT    = partT==i;             % Time points for this partition 
            idxQ    = partQ==i;             % Regressors for this partition 
            numFilt = size(xX.K(i).X0,2);   % Number of filter variables for this run 
            [Sw_hat(:,:,i),shrink(i)]=rsa.stat.covdiag(res(idxT,:),sum(idxT)-sum(idxQ)-numFilt-1);                    %%% regularize Sw_hat through optimal shrinkage
            [V,L]=eig(Sw_hat(:,:,i));       % This is overall faster and numerical more stable than Sw_hat.^-1/2
            l=diag(L);
            sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
            u_hat(idxQ,:)=beta_hat(idxQ,:)*sq;
        end;
        shrink=mean(shrink);
        
    case 'overall'
        [Sw_hat,shrink]=rsa.stat.covdiag(res,SPM.xX.erdf); %%% regularize Sw_hat through optimal shrinkage
        [V,L]=eig(Sw_hat);
        l=diag(L);
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        u_hat=beta_hat*sq;
end;
if (nargout>1)
    resMS=sum(res.^2)./SPM.xX.erdf;
end; 
if (nargout>2)
    Sw_raw=res'*res./SPM.xX.erdf;
end; 
