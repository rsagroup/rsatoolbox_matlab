% This fucntion masks any MNE source solution using a set of anatomically
% defined ROIs. The areas outside these ROIs are set to zeros. Note, this script
% also merges serveral ROIs into a bigger ROI. 
%
% Li Su 2-2012
% Note: This function assumes that all ROIs you give as input are for the same
% hemisphere. IZ

function maskedBrain = masking(brain, indexMasks, RoIs)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

allRegions = fieldnames(indexMasks); 
nMaskRegions = size(allRegions,1);
number_of_vertex = size(brain,1);

allRegions_lh = []; % changed these to arrays from cells to avoid any incosistencies in size % 03/12 IZ
allRegions_rh = [];

for RoI=1:numel(RoIs) %left

    thisMask = allRegions{RoIs(RoI)};
%     twFiledName = fieldnames(indexMasks.(thisMask));
    maskIndices = indexMasks.(thisMask).maskIndices;
    chi = indexMasks.(thisMask).chirality;
    maskIndices = maskIndices(maskIndices < number_of_vertex);
    
    if strcmp(chi, 'L')        
        allRegions_lh = [allRegions_lh maskIndices]; % changed to row vector IZ
    else
        allRegions_rh = [allRegions_rh maskIndices];
    end
    
end

allRegions_lh = sort(allRegions_lh);
allRegions_rh = sort(allRegions_rh);
wholeBrain = 1:number_of_vertex;

if strcmp(chi, 'L')     %% ?? which chi? the loopr above will have last chi stored only IZ 03/12
    zeroRegions = setdiff(wholeBrain',allRegions_lh);
else
    zeroRegions = setdiff(wholeBrain',allRegions_rh);    
end

maskedBrain = brain;
maskedBrain(zeroRegions,:) = 0;
