% Old: [RDM_fdtFolded_ltv, cv2RDM_fdt_sq] = fisherDiscrRDM(Xa, Ya, Xb, Yb, condPredIs, RDMmask)
function [RDM_fdcFolded_ltv, cv2RDM_fdc_sq] = fisherDiscrTRDM(Xa, Ya, Xb, Yb, condPredIs, RDMmask)
%
% USAGE
%       [RDM_fdcFolded, sdRDM_fdt] = fisherDiscrTRDM(Xa, Ya, Xb, Yb[, condPredIs, RDMmask)
%
%
%
% INPUT:
%   Xa : design matrix of first fold. nVolumes by nPreds
%   Ya : time series of first fold. nVolumes by nVoxels
%   Xb : design matrix of second fold. nVolumes by nPreds
%   Yb : time-series of second fold. 
%   condPredIs: vector [1 :  nConditions] indexing condition related elements of the OLS. 
% 
%
% computes a lower-triangular vector of Linear Discriminat Contrast
% (RDM_fdcFolded_ltv) and a folded symmetric LDC RDM (cv2RDM_fdc_sq).
% Inputs to this function are the design matrices and activation matrices
% for two independent datasets (a and b). CondPredIs is the indices for
% which the discriminabilities will be computed (usually 1:number of
% conditions). RDMmask constraints the analysis to the ones that are marked
% 1 only.
% LDC values are based on error covariance matrix. The shrinkage method
% introduced by Ledolt and Wolf (2013) is used to estimate the covariance
% matrix.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

% Modified August 2018 - Ian Charest - University of Birmingham.
%
% The function now returns the cross-validated mahalanobis distance. (LDC);
% the function now allows different number of runs (for example, splitting
% odd and even runs when total number of runs is odd.) This is achieved by creating 
% two separate contrast matrices C1 and C2, which match Xa and Xb's number 
% of columns respectively.

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

% preparations
nPred1=size(Xa,2);
nPred2=size(Xb,2);
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
C1=zeros(nPred1,nSelectedCondPairs); % nPred by nSelectedCondPairs
C1(condPredIs,:)=CC;

% IC: added a second contrast matrix. 
C2=zeros(nPred2,nSelectedCondPairs); % nPred by nSelectedCondPairs
C2(condPredIs,:)=CC;


%% fit Fisher dicriminant to set A, perform t test on set-B, and vice versa
% fishAtestB

[cs_ab]=fishAtestB_optShrinkageCov_C(Xa,Ya,Xb,Yb,C1,C2);
[cs_ba]=fishAtestB_optShrinkageCov_C(Xb,Yb,Xa,Ya,C2,C1);

% old
%[ps,ts_ab,c_ab]=fishAtestB_optShrinkageCov_C(Xa,Ya,Xb,Yb,C);
%[ps,ts_ba,c_ba]=fishAtestB_optShrinkageCov_C(Xb,Yb,Xa,Ya,C);


%% assemble the 2-fold cross-validation RDM (cv2RDM)
% results from the 2 folds in separate entries symmetric about the diagonal 
cv2RDM_fdc_sq=nan(nCond,nCond);
cv2RDM_fdc_sq(logical(RDMmask_ltsq))=cs_ab;
cv2RDM_fdc_sq=cv2RDM_fdc_sq';
cv2RDM_fdc_sq(logical(RDMmask_ltsq))=cs_ba;

%% fold the cv2RDM
% RDM_fdtFolded_ltv=foldCrossVal2RDM(cv2RDM_fdt_sq);

RDM_fdcFolded_ltv = mean(cat(2,cs_ab,cs_ba),2);


end%function


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunction: fishAtestB_optShrinkageCov_C
function [cs,invSa]=fishAtestB_optShrinkageCov_C(Xa,Ya,Xb,Yb,C1,C2,invSa)
% old: function [ps,ts,invSa]=fishAtestB_optShrinkageCov_C(Xa,Ya,Xb,Yb,C,invSa)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

numVox = size(Yb,2);

% fit model to data set A    
eBa=inv(Xa'*Xa)*Xa'*Ya; % nPredictorsA by nVox
eEa=Ya-Xa*eBa;          % compute residuals on set A (nTimePoints x nVoxels).

%% determine Fisher discriminant using data set A
if ~exist('invSa','var')||isempty(invSa)
    invSa=inv(covdiag(eEa)); % nVox by nVox covariance matrix of the error   
end
was=C1'*eBa*invSa; % noise normalised fisher dimension row vectors w (nComparisons by nVoxels) for data set A

%% test on data set B projected onto the set-A Fisher discriminant    
invXTXb=inv(Xb'*Xb); % (nPredictorsB x nPredictorsB)

yb_was=Yb*was'; % project time-series Yb (nTimePoints x nVoxels) onto wa (nComparisons by nVoxels) ==> (nTimePoints by nComparisons)
ebb_was=invXTXb*(Xb'*yb_was); % project invXTXb the invert design matrix b onto the fitted noise normalised Yb (nPredictorsB by nComparisons)

% I dropped this as we now use a contrast matrix from set b. this allows set a and set b having different number of runs.
%m1 = size(ebb_was,1);
%m2 = size(C,1);
%C_new = C(1:min([m1,m2]),:);% I CONSTRAINED C HERE.
%ctb_was2=diag(C_new'*ebb_was);                      % nContrasts by 1

cs=diag(C2'*ebb_was)/numVox;      
% C2'*ebb_was gives a (nComparisons x nComparisons) and its diagonal is the vector form of the LDC RDM.
% note that we now normalise this by the number of voxels.

%IC: I dropped the tRDMs and ps. Not quite sure how these pseudo-t-values (or ps) can be averaged anyways.
%se_ctb_was2=sqrt(esb_was.*diag(C_new'*invXTXb*C_new));  % 1 by nContrasts I CHANGED THIS FROM C TO C_NEW
%ts=ctb_was2./se_ctb_was2;                    % nContrasts by 1
%ps=tcdf(-ts,nDFb);                  % compute ps

end%function


%% utility anonymous functions below.

function sigma=covdiag(x)
% function [sigma]=covdiag(x)
%
% Regularises the estimate of the covariance matrix accroding to the 
% optimal shrikage method as outline in Ledoit& Wolf (2005). 
%
% Shrinks towards diagonal matrix
% INPUT: 
%       x     (t*n): t observations on p random variables
% OUTPUT: 
%       sigma (n*n): invertible covariance matrix estimator
 
% de-mean returns
[t,n]=size(x);
meanx=mean(x);
x=x-meanx(ones(t,1),:);
 
% compute sample covariance matrix
sample=(1/t).*(x'*x);
 
% compute prior
prior=diag(diag(sample));
 
% compute shrinkage parameters
d=1/n*norm(sample-prior,'fro')^2;
y=x.^2;
r2=1/n/t^2*sum(sum(y'*y))-1/n/t*sum(sum(sample.^2));
 
% compute the estimator
shrinkage=max(0,min(1,r2/d));
sigma=shrinkage*prior+(1-shrinkage)*sample;
end


function RDMs=squareRDMs(RDMs_ltv)
% converts set of row-vector RDMs_ltv to square form (despite being
% rows, RDMs are stacked along the 3rd dimension, just as square RDMs
% would be. this avoids ambiguity when the RDMs_ltv is square and could
% be either a single RDM or a number of vectorized RDMs.)
% RDMs may be bare or wrapped with meta-data in a struct object. they
% will be returned in the same format as passed.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council


if isstruct(RDMs_ltv)
    % wrapped
    RDMs_ltv_struct=RDMs_ltv;
    RDMs_ltv=unwrapRDMs(RDMs_ltv_struct);
    
    nRDMs=size(RDMs_ltv,3);
    RDMs=[];
    for RDMI=1:nRDMs
        RDMs=cat(3,RDMs,squareRDM(RDMs_ltv(:,:,RDMI)));
    end
    
    RDMs=wrapRDMs(RDMs,RDMs_ltv_struct);
else
    % bare
    nRDMs=size(RDMs_ltv,3);
    RDMs=[];
    for RDMI=1:nRDMs
        RDMs=cat(3,RDMs,squareform(vectorizeRDM(RDMs_ltv(:,:,RDMI))));
    end
end

end%function


function RDMs_struct=wrapRDMs(RDMs,refRDMs_struct)
% wraps similarity matrices RDMs (in square or upper triangle form)
% into a structured array with meta data copied from refRDMs_struct
% (which needs to have the same number of RDMs).(if they are already
% wrapped then the wrapping (metadata) is replaced by that of
% refRDMs_struct.
% generally in cases where one wants to wrap a number of dissimilarity
% vectors or matrices, they can define them as different 'RDM' fields of an
% input structure. In cases where the user wants to define names or specify
% colours for the dissimilarity vector or matrices they neeed to specify
% this in the refRDMs_struct, otherwise generic names and values would be
% selected and returned in the 'name' and 'color' fields.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

if isstruct(RDMs)
    % wrapped already, but replace the wrapping
    nRDMs=numel(RDMs);
else
    nRDMs=size(RDMs,3);
end    

if ~exist('refRDMs_struct','var')
    for RDMI=1:nRDMs
        refRDMs_struct(RDMI).name='[unnamed RDM]';
        refRDMs_struct(RDMI).color=[0 0 0];
    end
end

RDMs_struct=refRDMs_struct;
if isstruct(RDMs)
    % wrapped already, but replace the wrapping
    for RDMI=1:nRDMs
        RDMs_struct(RDMI).RDM=RDMs(RDMI).RDM;
    end
else
    % RDMs need wrapping
    for RDMI=1:nRDMs
        RDMs_struct(RDMI).RDM=RDMs(:,:,RDMI);
    end
end

end%function

function RDMs_utv=vectorizeRDMs(RDMs)
% converts set of RDMs (stacked along the 3rd dimension)
% to lower-triangular form (set of row vectors)
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

if isstruct(RDMs)
    % wrapped
    RDMs_struct=RDMs;
    RDMs=unwrapRDMs(RDMs_struct);
    
    nRDMs=size(RDMs,3);
    RDMs_utv=[];
    for RDMI=1:nRDMs
        RDMs_utv=cat(3,RDMs_utv,vectorizeRDM(RDMs(:,:,RDMI)));
    end
    
    RDMs_utv=wrapRDMs(RDMs_utv,RDMs_struct);
else
    % bare
    nRDMs=size(RDMs,3);
    RDMs_utv=[];
    for RDMI=1:nRDMs
        RDMs_utv=cat(3,RDMs_utv,vectorizeRDM(RDMs(:,:,RDMI)));
    end
end

end%function

function RDM=vectorizeRDM(RDM)
% converts RDM to lower-triangular form (row vector) (or leaves it in that
% form) in cases where the input(RDM) is a vector, the output would be the
% same as the input (row or column vector). if the input is a 2D matrix
% (e.g. an RDM) the output would be a (row) vector of the lower-triangular
% entries.
% this function does not accept a wrapped RDM as its input.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council


if size(RDM,1)==size(RDM,2)
    RDM(logical(eye(size(RDM))))=0; % fix diagonal: zero by definition
    RDM=squareform(RDM);
end

end%function



function [RDMs,nRDMs]=unwrapRDMs(RDMs_struct)
% unwraps dissimilarity matrices of a structured array (wrapped RDMs) with
% meta data by extracting the dissimilarity matrices (in square or lower
% triangle form) and lining them up along the third dimension. (if they are
% already in that format they are handed back unchanged.) it must be noted
% that different dissimilarity sets (i.e. lower-triangular or squareform
% RDMs) must be different fields of the input structure. The only
% requirement is that the dissimilarities should be stored in a field
% called 'RDM'.
%
% nRDMs is the number of wrapped RDMs in the input. in case where all the
% input RDMs are in lower-triangular or square form, the output would be in
% the same format as well (RDMs will be stacked along a third dimension).
% if the input RDMs are of mixed type, they're all made square. Also, if
% RDMs in the input are of different sizes, they'll all be returned as
% squares with nans outside. 
%
% Copyright (C) 2010 Medical Research Council

if isstruct(RDMs_struct)
    % in struct form
    nRDMs=size(RDMs_struct,2);
	
	mixedTypes = false;
	
	for RDMi = 1:nRDMs
		
		thisRDM = RDMs_struct(RDMi).RDM;
		
		% What type of RDM is this?
		if max(size(thisRDM)) == numel(thisRDM)
			% Then it's in ltv form
			thisType = 'ltv';
			thisSize = (1+sqrt(1-4*numel(thisRDM)))/2;
		else
			thisType = 'sq';
			thisSize = size(thisRDM,1);
		end%if:utv
		
		if RDMi == 1
			% What's the first type?
			firstType = thisType;
			firstSize = thisSize;
			maxSize = thisSize;
		elseif ~strcmp(thisType, firstType) || thisSize ~= firstSize;
			% Is each type the same as the first?
			mixedTypes = true;
			if thisSize ~= firstSize;
				maxSize = max(maxSize, thisSize);
			end%if:not the same size
		end%if:RDMi==1
	
		if ~mixedTypes
			switch thisType
				case 'ltv'
					RDMs(1,:,RDMi) = double(thisRDM);
				case 'sq'
					RDMs(:,:,RDMi) = double(thisRDM);
			end%switch:thisType
		end%if:~mixedTypes
		
	end%for:RDMi
	
	% If previous loop was broken...
	
	if mixedTypes
		clear RDMs;
		% All must be square
		RDMs = nan(maxSize, maxSize, nRDMs);
		for RDMi = 1:nRDMs
			thisRDM = RDMs_struct(RDMi).RDM;
			% What type of RDM is this?
			if max(size(thisRDM)) == numel(thisRDM)
				% Then it's in ltv form
				thisRDM = squareform(thisRDM);
			end%if:utv
			RDMs(1:size(thisRDM,1), 1:size(thisRDM,2), RDMi) = double(thisRDM);
		end%for:RDMi
	end%if:mixedTypes
	
else
    % bare already
    RDMs=RDMs_struct;
    nRDMs=size(RDMs,3);
end

end%function

