function [u_hat,resMS,Sw_hat,beta_hat,shrinkage,trRR]=noiseNormalizeBeta(Y,SPM,varargin)
% function [u_hat,Sw_hat,resMS,beta_hat]=rsa_noiseNormalizeBeta(Y,SPM,varargin)
% Estimates beta coefficiencts beta_hat and residuals from raw time series Y
% Estimates the true activity patterns u_hat by applying noise normalization to beta_hat
% INPUT:
%    Y        raw timeseries, T by P
%    SPM:     SPM structure
% OPTIONS:
%   'normmode':   'overall': Does the multivariate noise normalisation overall (default)
%                 'runwise': Does the multivariate noise normalisation by run
%   'shrinkage':  Shrinkage coefficient. 
%                 0: No regularisation 
%                 1: Using only the diagonal - i.e. univariate noise normalisation 
%                 By default the shrinkage coeffcient is determined using
%                 the Ledoit-Wolf method. 
% OUTPUT:
%    u_hat:       estimated true activity patterns (beta_hat after multivariate noise normalization),
%    resMS:       residual mean-square - diagonal of the Var-cov matrix of the average beta weight, 1 by P  
%    Sw_hat:      overall voxel error variance-covariance matrix (PxP), before regularisation 
%    beta_hat:    estimated raw regression coefficients, K*R by P
%    shrinkage:   applied shrinkage facto 
%    trRR:        Trace of the squared residual spatial correlation matrix
%                 This number provides an estimate for the residual spatial
%                 correlation. For variance-estimation, the effective
%                 number of voxels are effVox = numVox^2/trRR 
% Alexander Walther, Joern Diedrichsen
% joern.diedrichsen@googlemail.com
% 10/2018: Updated scaling of prehwitened-betas, to take into account the
% variance of betas (bCov).
Opt.normmode = 'overall';  % Either runwise or overall
Opt.shrinkage = []; 
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

%%% redo the first-level GLM using matlab functions 
KWY=spm_filter(xX.K,xX.W*Y);                               %%% filter out low-frequence trends in Y
beta_hat=xX.pKX*KWY;                                       %%% ordinary least squares estimate of beta_hat = inv(X'*X)*X'*Y
res=spm_sp('r',xX.xKXs,KWY);                               %%% residuals: res  = Y - X*beta
clear KWY XZ                                               %%% clear to save memory

switch (Opt.normmode)
    case 'runwise'              % do run-wise noise normalization
        u_hat   = zeros(size(beta_hat));
        shrink=zeros(Nrun,1);
        for i=1:Nrun
            idxT    = partT==i;             % Time points for this partition 
            idxQ    = partQ==i;             % Regressors for this partition 
            numFilt = size(xX.K(i).X0,2);   % Number of filter variables for this run 
            
            % in the scaling of the noise, take into account mean beta-variance 
            [Sw_reg(:,:,i),shrinkage(i),Sw_hat(:,:,i)]=rsa.stat.covdiag(res(idxT,:),SPM.xX.trRV/(Nrun*mean(diag(SPM.xX.Bcov))),'shrinkage',Opt.shrinkage);                    %%% regularize Sw_hat through optimal shrinkage
            % Calculating sq over the eigenvalues is numerically more
            % stable than sq = Sw_reg^-1/2 
            [V,L]=eig(Sw_reg(:,:,i));   
            l=diag(L);
            sq(:,:,i) = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
            % Postmultiply by the inverse square root of the estimated matrix 
            u_hat(idxQ,:)=beta_hat(idxQ,:)*sq(:,:,i);
        end;
        shrinkage=mean(shrinkage);
        Sw_hat = mean(Sw_hat,3); 
        Sw_reg = mean(Sw_reg,3); 
    case 'overall'              %%% do overall noise normalization
        % in the scaling of the noise, take into account mean beta-variance 
        [Sw_reg,shrinkage,Sw_hat]=rsa.stat.covdiag(res,SPM.xX.trRV/mean(diag(SPM.xX.Bcov)),'shrinkage',Opt.shrinkage);   %%% regularize Sw_hat through optimal shrinkage

        % Postmultiply by the inverse square root of the estimated matrix 
        [V,L]=eig(Sw_reg);  % in the scaling of the noise, take into account mean beta-variance 
        l=diag(L);
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        u_hat=beta_hat*sq;
end;

% Return the diagonal of Sw_hat - also weighted by the mean variance of beta-hat 
if (nargout>1)
    resMS=sum(res.^2)./SPM.xX.trRV*mean(diag(SPM.xX.Bcov));
end; 

% return the trace of the residual covariance matrix
% If prewhiting was perfect, then the trace(sig*sig) = P and the effective number of voxels
% would be P. However, because we need to regularise our estimate, the residual 
% spatial covariance matrix is somwhat correlated. For variance estimation on the distances, 
% the effective voxel number would be numVox^2/trRR
if (nargout>5)
    % Caluclate the predicted redidual covariance 
    sq     = mean(sq,3); 
    V_hat  = sq'*Sw_hat*sq; 
    % Calculate the correlation matrix from this 
    Vsigma = sqrt(diag(V_hat));     
    R_hat = bsxfun(@rdivide,V_hat,Vsigma); 
    R_hat = bsxfun(@rdivide,R_hat,Vsigma'); % R = C ./ sigma*sigma';
    % get the trace 
    trRR  =  sum(sum(R_hat.*R_hat)); % traceRR 
end; 
