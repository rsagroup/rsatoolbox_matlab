function [RDM_fdtFolded_ltv, cv2RDM_fdt_sq] = fisherDiscrTRDM(Xa, Ya, Xb, Yb, condPredIs, RDMmask)

% USAGE
%       [RDM_fdtFolded, sdRDM_fdt] = fisherDiscrTRDM(Xa, Ya, Xb, Yb[, condPredIs, RDMmask)
%
%
% computes a lower-triangular vector of Linear Discriminat t-statistics
% (RDM_fdtFolded_ltv) and a folded symmetric LDt RDM (cv2RDM_fdt_sq).
% Inputs to this function are the design matrices and activation matrices
% for two independent datasets (a and b). CondPredIs is the indices for
% which the discriminabilities will be computed (usually 1:number of
% conditions). RDMmask constraints the analysis to the ones that are marked
% 1 only.
% LDt values are based on error covariance matrix. The shrinkage method
% introduced by Ledolt and Wolf (2013) is used to estimate the covariance
% matrix.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% preparations
nPred=size(Xa,2);
if ~exist('condPredIs','var'), condPredIs=1:nPred; end
nCond=numel(condPredIs);
if ~exist('RDMmask','var'), RDMmask=true(nCond,nCond); end
RDMmask_ltsq=tril(squareRDMs(RDMmask),-1); % square RDM mask restricted to lower triangular region


%% construct contrasts matrix C
[i,j]=find(RDMmask_ltsq);
nSelectedCondPairs=length(i);
CC=zeros(nCond,nSelectedCondPairs); % nCond by nSelectedCondPairs
ii=sub2ind(size(CC),i,(1:length(i))');
jj=sub2ind(size(CC),j,(1:length(j))');
CC(ii)=-1;
CC(jj)=1;
C=zeros(nPred,nSelectedCondPairs); % nPred by nSelectedCondPairs
C(condPredIs,:)=CC;


%% fit Fisher dicriminant to set A, perform t test on set-B, and vice versa
% fishAtestB

[ps,ts_ab]=fishAtestB_optShrinkageCov_C(Xa,Ya,Xb,Yb,C);
[ps,ts_ba]=fishAtestB_optShrinkageCov_C(Xb,Yb,Xa,Ya,C);

%% assemble the 2-fold cross-validation RDM (cv2RDM)
% results from the 2 folds in separate entries symmetric about the diagonal 
cv2RDM_fdt_sq=nan(nCond,nCond);
cv2RDM_fdt_sq(logical(RDMmask_ltsq))=ts_ab;
cv2RDM_fdt_sq2=cv2RDM_fdt_sq';
cv2RDM_fdt_sq2(logical(RDMmask_ltsq))=ts_ba;

%% fold the cv2RDM
% RDM_fdtFolded_ltv=foldCrossVal2RDM(cv2RDM_fdt_sq);
cv2RDM_fdt_sq(find(isnan(cv2RDM_fdt_sq)))=0;
RDM_fdtFolded_ltv = squareform(cv2RDM_fdt_sq);

end%function


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunction: fishAtestB_optShrinkageCov_C
function [ps,ts,invSa]=fishAtestB_optShrinkageCov_C(Xa,Ya,Xb,Yb,C,invSa)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% fit model to data set A    
eBa=inv(Xa'*Xa)*Xa'*Ya; % nCond by nVox
eEa=Ya-Xa*eBa;

%% determine Fisher discriminant using data set A
if ~exist('invSa','var')||isempty(invSa)
    invSa=inv(covdiag(eEa)); % nVox by nVox
end
was=C'*eBa*invSa; % fisher dimension row vectors w (nContrasts by nVox) for data set A

%% t test on data set B projected onto the set-A Fisher discriminant    
invXTXb=inv(Xb'*Xb);

yb_was=Yb*was'; % project TCs Yb onto wa (nTimePoints by nContrasts)
ebb_was=invXTXb*(Xb'*yb_was);             % nPredictors by nContrasts
eeb_was=yb_was-Xb*ebb_was;              % nTimePoints by nContrasts
nDFb=size(yb_was,1)-size(Xb,2);         % n degrees of freedom
esb_was=diag(eeb_was'*eeb_was)/nDFb;            % error variance on wa (scalar for each contrast: nContrasts by 1)

m1 = size(ebb_was,1);
m2 = size(C,1);
C_new = C(1:min([m1,m2]),:);% I CONSTRAINED C HERE.
ctb_was2=diag(C_new'*ebb_was);                      % nContrasts by 1

se_ctb_was2=sqrt(esb_was.*diag(C_new'*invXTXb*C_new));  % 1 by nContrasts I CHANGED THIS FROM C TO C_NEW
ts=ctb_was2./se_ctb_was2;                    % nContrasts by 1

ps=tcdf(-ts,nDFb);                  % compute ps

end%function
