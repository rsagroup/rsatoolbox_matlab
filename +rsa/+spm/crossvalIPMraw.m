function G = crossvalIPMraw(Y,SPM,conditionVec,varargin)
% function dist=rsa.spm.crossvalIPMraw(Y,SPM,conditionVec,varargin);
% First, gets the regression coefficent from the SPM, and prewhitens them.
% Prewhiten can be controlled using different methods (run-wise or overall).
% The same code is used in rsa.spm.noiseNormalizeBeta.
% Secondly, it calcualtes inner product (second moment) matrix with a 
% leave-one out crossvalidation,
% It uses optimal combination of the beta-coeeficients in the part that
% averages across partitions. By default, the different partitions are
% assumed to be the different imaging runs.
% Note that in calculating the inner product, the structure of the first-level 
% design matrix is optimally taken into account. The command: 
%  IPM  = rsa.spm.crossvalIPMraw(Y,SPM,condition); 
% is therefore equivalent to:  
%  beta = rsa.spm.noiseNormalizeBeta(Y,SPM); 
%  IPM  = rsa.crossvalIPM(beta,partition,condition,SPM.xX.xKXs.X); 
% INPUT:
%    Y             raw timeseries, T by P
%    SPM:          SPM structure
%    conditionVec: Vector that indicates conditions for each column of design matrix 
%                  set to 0 for nuisance regressors
% OUTPUT:
%    IPM:          numCondxnumCond inner product matrix
% 
% OPTIONs:
%   'normmode':    'runwise': Does the multivariate noise normalisation by run
%                  
%                  'overall': Does the multivariate noise normalisation overall
%   'normmethod':  'multivariate': The is the default using ledoit-wolf reg.
%                  'univariate': Performing univariate noise normalisation (t-values)
%                  'none':    No noise normalisation
% 
%   'partition':    In case you want to use specific partitioning, pass a structure
%                   with fields; 'partT' and 'partN'
% 
% (c) 2015 Joern Diedrichsen, Alex Walther
Opt.normmode = 'runwise';  % Either runwise or overall
Opt.normmethod = 'multivariate';  % Either runwise or overall
Opt.partition = [];
Opt = rsa.getUserOptions(varargin,Opt,...
    {'normmode','normmethod','partition'});

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

%%% Get partions: For each run (1:K), find the time points (T) and regressors (K+Q) that belong to the run
partT = nan(T,1);
partN = nan(numReg,1);
numPart=length(SPM.Sess);                                     %%% number of runs
for i=1:numPart
    partT(SPM.Sess(i).row,1)=i;
    partN(SPM.Sess(i).col,1)=i;
    partN(SPM.xX.iB(i),1)=i;                                %%% Add intercepts
end;

if ~isempty(Opt.partition); 
    partT = Opt.partition.partT; 
    partN = Opt.partition.partN;
    numPT = length(unique(partT));
    numPN = length(unique(partN));
    if numPT~=numPN
        error('Number of partition is not equal!');
    else
        numPart = numPT;
    end;    
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

[G] = getCVIPM(nonInterest,partN,partT,Z,X,KWY,numVox,[1:numPart]);

end

% Local function to get cross-validated distance
function IPM = getCVIPM(nonInterest,partN,partT,Z,X,KWY,numVox,partition)
    nPartition = length(partition);
    A = zeros(length(nonInterest(partN==partN(1))),numVox,nPartition);
    for p=1:nPartition
        % Get the betas from the test run
        indxN = partN==partition(p);
        indxT = partT==partition(p);
        Za = Z(indxN,:);
        Za = Za(:,any(Za,1));
        Xa = X(indxT,indxN);
        Ma  = Xa*Za;
        A(:,:,p)     = (Ma'*Ma)\(Ma'*KWY(indxT,:));
        
        % Get the betas based on the other runs
        indxN = partN~=partition(p);
        indxT = partT~=partition(p);
        Zb    = Z(indxN,:);
        Zb    = Zb(:,any(Zb,1));
        Xb    = X(indxT,indxN);
        Mb    = Xb*Zb;
        B     = (Mb'*Mb)\Mb'*KWY(indxT,:);
        
        % Pick condition of interest
        interest = ~nonInterest(partN==partition(p));
        IPM(:,:,p) = A(interest,:,p)*B(interest,:)'/numVox;        
    end;
    IPM = mean(IPM,3);
    IPM = 0.5*(IPM+IPM');
end