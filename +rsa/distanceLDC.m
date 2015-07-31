function d=distanceLDC(Y,partition,conditionVec)
% function d=distanceLDC(Y,partition,conditionVec);
% Calculates the crossvalidated squared Euclidian distances 
% Trial on 
% INPUT:
%  Y           : Noise-normalized activation patterns, a KR x P matrix
%  partition   : KR x 1 integer value that indicates the partition for crossvalidation (typically run number)
%  conditionVec: KR x 1 vector of conditions, zeros will be ignored 
% 
% OUTPUT:
%   d         : Average distance values for each contrast across folds
%               (1xC), this measure is normalised to the number of voxels 
% Alexander Walther, Joern Diedrichsen 
% 2/2015 

import rsa.util.*; 

[N,numVox]   = size(Y); 
part    = unique(partition)';
numPart = numel(part);
numCond   = max(conditionVec); 

A = zeros(numCond,numVox,numPart); 
X = indicatorMatrix('identity_p',conditionVec); 
C = indicatorMatrix('allpairs',[1:numCond]); 

% Estimate condition means within each run 
for i=1:numPart 
    Xa = X(partition==part(i),:);
    Ya = Y(partition==part(i),:);
    Xb = X(partition~=part(i),:);
    Yb = Y(partition~=part(i),:);
    A(:,:,i) = pinv(Xa)*Ya;
    B        = pinv(Xb)*Yb; 
    d(i,:)   = sum((C*A(:,:,i)).*(C*B),2)'/numVox;      % Note that this is normalised to the number of voxels 
end; 
d = sum(d)./numPart; 

