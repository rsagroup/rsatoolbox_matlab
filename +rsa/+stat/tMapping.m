function [tVecs,eB,eE]=tMapping(Y,X,C)

% USAGE
%       [tVecs,eB,eE]=tMapping(Y,X,C)
%
% FUNCTION
%       on the basis of the time-by-space data matrix Y, this funtion
%       computes map of t statistics characterizing the effects defined by
%       contrast column vectors C referring to design matrix X.
%
%       if C is a matrix with more than one column, a t map is computed for
%       each column of C.
%
%       the t maps are returned as space vectors tVecs: a space-by-contrast
%       matrix.
%
% ARGUMENTS
% Y             the 3D+time fMRI data set as a time-by-space
%               matrix.
%
% X             the time-by-predictor design matrix corresponding to data
%               set Y.
%
% C             the contrast column vector referring to the design matrix X
%               (length(C) must be size(X,2)) which defined the contrast
%               the fMRI volume is to t-mapped for. if C is a matrix, then
%               a t map is computed for each column vector (defining a
%               contrast) and the returned tVecs matrix will have a
%               second dimension (voxel by contrast).

[nPreds,nContrasts]=size(C);
[nTime,nVox]=size(Y);

tVecs=zeros(nVox,nContrasts);

XtXi=inv(X'*X);                 % model-by-model

eB=XtXi*X'*Y;                   % model-by-space: (e)stimate of B (beta matrix) minimizing sum of squared deviations

nDF=nTime-nPreds;        % scalar: (n)umber of (d)egrees of (f)reedom

% do the following with less memory...
% eE=Y-X*eB;                      % time-by-space: (e)stimated of E (errors matrix)
% evarE=sum(eE.^2,1)/nDF;         % 1-by-space: vertical integration resulting in row vector containing one (e)stimate of the (var)iance of the errors (E) for each voxel

eE=X*eB;                        % time-by-model
evarE=zeros(1,size(Y,2));
for timeI=1:size(Y,1) % for each time point...
    eE(timeI,:)=Y(timeI,:)-eE(timeI,:);
    evarE=evarE+eE(timeI,:).^2;         % 1-by-space: vertical integration resulting in row vector containing one (e)stimate of the (var)iance of the errors (E) for each voxel
end
evarE=evarE/nDF;

for contrastI=1:nContrasts
    eeffSTD=sqrt(C(:,contrastI)'*XtXi*C(:,contrastI)*evarE);        % 1-by-space: (e)stimated (st)andard (d)eviation of the estimated (eff)ect defined by contrast vector c
    tVecs(:,contrastI)=C(:,contrastI)'*eB./eeffSTD;                                  % 1-by-space: row vector of t effects
end


