function d=distanceLDCraw(Y,SPM,conditionVec)
% function d=rsa.spm.distanceLDCraw(Y,SPM,conditionVec); 
% First, gets the regression coefficent using the same methods as 
% rsa.spm.noiseNormalizeBeta (prewhitening within each imaging run
% seperately) 
% Secondly, it calcualtes LDC stances with a leave-one out crossvalidation,
% It uses optimal combination of the beta-coeeficients in the part that
% averages across partitions. By default, the different partitions are
% assumed to be the different imaging runs. 
% INPUT: 
%    Y            raw timeseries, T by P 
%    SPM:         SPM structure
%    conditionVec: a SPM.xX.iC-long ConditionVector, that indicates
%                 conditions, set 0 for nuisance regressors 
% OUTPUT: 
%    d:           numCond*(numCond-1)/2 distances between experimental
%                 conditions 
% Joern Diedrichsen, Alex Walther 

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
for i=1:numPart
    idxT=partT==i;
    Sw_hat(:,:,i)=covdiag(res(idxT,:));                    %%% regularize Sw_hat through optimal shrinkage
    KWY(idxT,:)=KWY(idxT,:)*Sw_hat(:,:,i)^(-1/2);   %%% multivariate noise normalization
end;

% Estimate condition means within each 
for i=1:numPart 
    p     = [1:numPart]; 
    indxN = partN==i;
    indxT = partT==i;
    Za = Z(indxN,:); 
    Za = Za(:,any(Za,1)); 
    Xa = X(indxT,indxN);
    Ma  = Xa*Za; 
    A     = (Ma'*Ma)\(Ma'*KWY(indxT,:));
    indxN = partN~=i;
    indxT = partT~=i;
    p(i)  = []; 
    Zb    = Z(indxN,:); 
    Zb    = Zb(:,any(Zb,1)); 
    Xb    = X(indxT,indxN);
    Mb    = Xb*Zb; 
    B     = (Mb'*Mb)\Mb'*KWY(indxT,:);
    d(i,:)= sum((C*A(1:numCond,:)).*(C*B(1:numCond,:)),2)'/numVox;      % Note that this is normalised to the number of voxels 
end; 
d = sum(d)./numPart; 

