function [u_hat,resMS,Sw_hat,beta_hat]=noiseNormalizeBeta(Y,SPM)
% function [u_hat,Sw_hat,resMS,beta_hat]=rsa_noiseNormalizeBeta(Y,SPM)
% % % Estimates beta coefficiencts beta_hat and residuals from raw time series Y 
% % % Estimates the true activity patterns u_hat by applying multivariate noise
% % % normalization to beta_hat
% % % INPUT: 
% % %    Y        raw timeseries, T by P 
% % %    SPM:     SPM structure
% % % OPTIONS:
% % %    condidx  indices of the experimental conditions
% % % OUTPUT: 
% % %    u_hat    estimated true activity patterns (beta_hat after multivariate noise normalization),
% % %    resMS    residual mean-square - diagonal of the Var-cov matrix, 1 by P 
% % %    Sw_hat   estimated voxel error variance-covariance matrix of each run, P by P
% % %    beta_hat estimated regression coefficients, K*R by P
% % % Alexander Walther, Joern Diedrichsen 
% % % joern.diedrichsen@googlemail.com
% % % 2/2015

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

resMS=sum(res.^2)./(T-Q); 

%%% do run-wise noise normalization  
for i=1:Nrun
    idxT=partT==i;
    idxQ=partQ==i;
    Sw_hat(:,:,i)=covdiag(res(idxT,:));                    %%% regularize Sw_hat through optimal shrinkage
    u_hat(idxQ,:)=beta_hat(idxQ,:)*Sw_hat(:,:,i)^(-1/2);   %%% multivariate noise normalization
end;