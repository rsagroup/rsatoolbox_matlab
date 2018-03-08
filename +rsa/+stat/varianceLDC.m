function V=varianceLDC(d,C,Sig,nPart,nVox); 
% function V=varianceLDC(d,C,Sig,nPart,nVox); 
% returns predicted variance-covariance of a set of LDC distances, which are 
% assumed to caluclated using rsa.distanceLDC. 
% INPUT: 
%   d/G:    d: (nDist x 1) or (1 x nDist) vector of expected distances, 
%               these are assumed to be in a squareform(pdist) order  
%           G: A nCond x nCond second momement matrix of the patterns 
%   C:      Contrast matrix that gives rise to the distances (nDist x
%           nCond). If empty, it assumes that all pairwise distances are caluclated 
%           this uses optimised code. 
%   Sig:    Variance-covariance matrix of noise in each fold of the
%           original coefficicents (nCond x nCond), or scalar, assuming Sig = eye(nCond)*sig
%           Noise is also assumed to be independent across voxels
%   nPart:  Number of independent partitions
%   nVox:   (effective) number of voxels 
% OUTPUT: 
%   V:   Variance-covariance matrix of the distances, a (DxD) matrix 
% (c) 2015 Joern Diedrichsen (joern.diedrichsen@googlemail.com)

% Don't check input sizes too avoid overhead 
nFolds        = nPart*(nPart-1);        % Number of pairs of partitions

% Infer the original G to get variance-covariance of the distances 
[m,n]=size(d); 
if (m==n) % Square matrix: Second moment 
    G=d; 
else 
    G  =  squareform(-0.5*d); % Vector: distances
end; 
nCond = size(G,1); 

% Get the within-fold variance-covariance matrix of distances 
if (isscalar(Sig)) 
    Sig = eye(nCond).*Sig;                    
end; 
if (isempty(C))
    dG = rsa.stat.pairCACt(G);                         % Predicted variance-covariance of the signals
    dSig = rsa.stat.pairCACt(Sig); 
else 
    dG = C*G*C'; 
    dSig = C*Sig*C'; 
end; 

% Now Calculate the variance-covatiance matrix of the distances 
V=(4*(dG.*dSig)/nPart + 2*(dSig.*dSig)/nFolds)/nVox;
