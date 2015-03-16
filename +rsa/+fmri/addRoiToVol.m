function vol=addRoiToVol(vol, roi, color, solid)

% USAGE:        vol=addRoiToVol(vol, roi[, color, solid])
%
% FUNCTION:     to mark the specified region of interest to the true-color
%               volume vol.
%
% vol           the volume as a stack of true-color slices: X by Y by 3 by Z.
%               the third dimension encodes the color component (red, green,
%               blue).
%               anatomically, the X axis points to the left, the Y axis to
%               the back of the brain, and the Z axis up.
%
% roi           a region of interest defined as a set of voxels. roi is a
%               matrix of size # voxels by 3, every row of which specifies
%               a voxel by its ONE-BASED coordinate triple (x,y,z).
%
% color         optional red-green-blue triple [r, g, b] specifying the color. values 
%               range from 0 to 1. defaults to cyan.
%
% solid         optional flag that determines whether the roi is marked
%               solidly (1, default) or only by its 3D contour (0) or by its
%               2D (within-slice) contour (2).
%
% GENERAL IDEA: If we already have a volume containing colour information
%               and we want to highlight an area (roi) in the volume (vol)
%               with a specific colour; then we indicate the new area by
%               the indices in "roi" and the new colour in "color".  A new
%               "vol" will be returned which has the area highlighted.  The
%               colour content of the other voxels (the ones outside roi)
%               remain unchanged in the new vol.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if ~exist('color','var')
    color=[0 1 1];
    %color=randomColor;
end

if ~exist('solid','var')
    solid=1;
end

if solid~=1
    sizeVol=[size(vol,1),size(vol,2),size(vol,4)];
    nSideNeighbors=zeros(sizeVol);

    roiVOL=zeros(sizeVol);        roiVOL(     sub2ind(sizeVol, roi(:,1), roi(:,2), roi(:,3) ) )=1;
    
    up=   [roi(:,1),roi(:,2),roi(:,3)+1];
    down= [roi(:,1),roi(:,2),roi(:,3)-1];
    back= [roi(:,1),roi(:,2)+1,roi(:,3)];
    forth=[roi(:,1),roi(:,2)-1,roi(:,3)];
    left= [roi(:,1)+1,roi(:,2),roi(:,3)];
    right=[roi(:,1)-1,roi(:,2),roi(:,3)];
   
    % exclude out-of-volume voxels
    outgrowths = up(:,1)<1 | up(:,2)<1 | up(:,3)<1 | up(:,1)>sizeVol(1) | up(:,2)>sizeVol(2) | up(:,3)>sizeVol(3);
    up(find(outgrowths),:)=[];

    outgrowths = down(:,1)<1 | down(:,2)<1 | down(:,3)<1 | down(:,1)>sizeVol(1) | down(:,2)>sizeVol(2) | down(:,3)>sizeVol(3);
    down(find(outgrowths),:)=[];

    outgrowths = back(:,1)<1 | back(:,2)<1 | back(:,3)<1 | back(:,1)>sizeVol(1) | back(:,2)>sizeVol(2) | back(:,3)>sizeVol(3);
    back(find(outgrowths),:)=[];
    
    outgrowths = forth(:,1)<1 | forth(:,2)<1 | forth(:,3)<1 | forth(:,1)>sizeVol(1) | forth(:,2)>sizeVol(2) | forth(:,3)>sizeVol(3);
    forth(find(outgrowths),:)=[];

    outgrowths = left(:,1)<1 | left(:,2)<1 | left(:,3)<1 | left(:,1)>sizeVol(1) | left(:,2)>sizeVol(2) | left(:,3)>sizeVol(3);
    left(find(outgrowths),:)=[];

    outgrowths = right(:,1)<1 | right(:,2)<1 | right(:,3)<1 | right(:,1)>sizeVol(1) | right(:,2)>sizeVol(2) | right(:,3)>sizeVol(3);
    right(find(outgrowths),:)=[];
    
    upVOL=zeros(sizeVol);    upVOL(      sub2ind(sizeVol, up(:,1), up(:,2), up(:,3) ) )=1;
    downVOL=zeros(sizeVol);  downVOL(    sub2ind(sizeVol, down(:,1), down(:,2), down(:,3) ) )=1;
    backVOL=zeros(sizeVol);  backVOL(    sub2ind(sizeVol, back(:,1), back(:,2), back(:,3) ) )=1;
    forthVOL=zeros(sizeVol); forthVOL(   sub2ind(sizeVol, forth(:,1), forth(:,2), forth(:,3) ) )=1;
    leftVOL=zeros(sizeVol);  leftVOL(    sub2ind(sizeVol, left(:,1), left(:,2), left(:,3) ) )=1;
    rightVOL=zeros(sizeVol); rightVOL(   sub2ind(sizeVol, right(:,1), right(:,2), right(:,3) ) )=1;
   
    if solid==0
        nSideNeighbors=upVOL+downVOL+backVOL+forthVOL+leftVOL+rightVOL;
        [ivolx,ivoly,ivolz]=ind2sub(sizeVol,find(roiVOL&(nSideNeighbors<6)));
    else % e.g. solid==2
        nSideNeighbors=backVOL+forthVOL+leftVOL+rightVOL;
        [ivolx,ivoly,ivolz]=ind2sub(sizeVol,find(roiVOL&(nSideNeighbors<4)));
    end
   
    roi=[ivolx,ivoly,ivolz];
end

if size(roi,1)==0;
    return;
end
    
for i=1:size(roi,1)
    vol(roi(i,1),roi(i,2),:,roi(i,3))=color;
end

% different meaning (nonsensical in this case)
% i=1:size(roi,1);
% vol(roi(i,1),roi(i,2),:,roi(i,3))=color;

end%function