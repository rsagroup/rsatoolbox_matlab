function pThreshold=FDRthreshold(pMap,q,binaryBrainMask,assumePositiveDependence)

% USAGE         pThreshold=FDRthreshold(pMap
%                                       [,q=0.05]
%                                       [,binaryBrainMask=ones(size(pMap))]
%                                       [,assumePositiveDependence=1])
%
% FUNCTION      to determine the p value, at which the 3D map of p values
%               pMap needs to thresholded in order for the average false
%               discovery rate (FDR) to be below q.
%
% ARGUMENTS 
% pMap          the 3D map of p values
%
% [q]           the average false discovery rate the researcher is willing
%               to accept. (the actual false discovery rate can be larger,
%               but its average is guaranteed to be smaller than q.)
%               this argument is optional and defaults to 0.05.
%
% [binaryBrainMask] a binary 3D map whose size must match pMap. only voxels
%               marked by a nonzero entry in binaryBrainMask are considered
%               in determining the threshold. specifying this mask is
%               highly recommended, because it reduces the expected number of
%               falsely positive voxels under the null hypothesis for any
%               given threshold and increases the sensitivity with which
%               functional regions are detected. however, this argument is
%               optional. if it is omitted the whole volume will be
%               considered.
%
% [assumePositiveDependence] optional argument specifying whether no
%               assumptions at all or very weak assumptions are to be made
%               about the joint distribution of the p values. if
%               assumePositiveDependence is 1 (default), it is assumed that
%               the p values are either INDEPENDENT or POSITIVELY DEPENDENT
%               (Benjamini and Yekutieli, 2001), i.e. the noise in the data
%               is gaussian with nonnegative correlation across voxels.
%               this assumption is usually reasonable for fMRI data. it
%               also appears reasonable to me for p values stemming from
%               local-pattern-effect mapping. if assumePositiveDependence
%               is 0, no assumptions are made about the joint distribution
%               of the p values. this comes at a cost in sensitivity, i.e.
%               the resulting thresholds will be smaller, marking fewer
%               voxels. by default, assumePositiveDependence is 1 and the
%               weak assumptions (as described) are made.
%
% RETURN VALUES
% pThreshold    the critical p value, at which pMap is to be thresholded in
%               order to ensure that the average false discovery rate does
%               not exceed q.
%
% AUTHOR        nikolaus kriegeskorte, 2005
%
% REFERENCE     Genovese CR, Lazar NA, and Nichols T. Thresholding of
%               statistical maps in functional neuroimaging using the false
%               discovery rate. NeuroImage 15, 870-878. 2002.

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% PREPARATIONS
if ~exist('assumePositiveDependence','var')
    assumePositiveDependence=1;
end

if ~exist('binaryBrainMask','var') || (exist('binaryBrainMask','var') && isempty(binaryBrainMask))
    binaryBrainMask=true(size(pMap));
end
binaryBrainMask=logical(binaryBrainMask);

if ~exist('q','var')
    q=0.05;
end


pMap_vec=pMap(binaryBrainMask);

if any(isnan(pMap_vec))
    disp('FDRthreshold: found NaNs and excluded them from the analysis.');
    pMap_vec=pMap_vec(~isnan(pMap_vec));
end


n=length(pMap_vec);

if assumePositiveDependence
    c=1;
else
    c=sum([1:n].^-1);
end


%% STEP 1: order the p values

pMap_vec_sorted=sort(pMap_vec);
pMap_vec_sorted = pMap_vec_sorted(:);

%% STEP 2: find the critical p value that yields an average FDR smaller than q
pThreshold = pMap_vec_sorted(max(find(pMap_vec_sorted<=(1:n)'/n*q/c)));
if ~isempty(pThreshold)
    nSig = numel(find(pMap_vec_sorted <= pThreshold));
else
    nSig = 0;
end

disp([num2str(nSig),' tests exceed the FDR threshold for q<',num2str(q),'.']);

if isempty(pThreshold)
    pThreshold = inf;
end

end%function
