function relRoi = sphericalRelativeRoi(rad_mm, voxelSize_mm);

% USAGE:        sphereRelRoi = sphericalRelativeRoi(rad)
%
% FUNCTION:     create a solid sphere of voxels of radius rad.
%
% ARGUMENTS:
% rad_mm        the sphere's radius in millimeters
%
% voxelSize_mm  a triple specifying the voxel size along the three
%               dimensions in the order x, y, z in millimeters.
%
% RETURN VALUES:
% sphereRelRoi  the voxel sphere as a roi matrix (nVox by 3),
%               i.e. a list of voxel locations.
%               each row [x, y, z] contains the RELATIVE
%               coordinates (center of the sphere at (0, 0, 0)).
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

rightHalf0=0:voxelSize_mm(1):rad_mm; leftHalf=fliplr(rightHalf0(2:size(rightHalf0,2)));
xrange=[leftHalf,rightHalf0];

rightHalf0=0:voxelSize_mm(2):rad_mm; leftHalf=fliplr(rightHalf0(2:size(rightHalf0,2)));
yrange=[leftHalf,rightHalf0];

rightHalf0=0:voxelSize_mm(3):rad_mm; leftHalf=fliplr(rightHalf0(2:size(rightHalf0,2)));
zrange=[leftHalf,rightHalf0];

[X_mm, Y_mm, Z_mm]=meshgrid(xrange, yrange, zrange);

voxelRadii=((X_mm.^2)+(Y_mm.^2)+(Z_mm.^2)).^(0.5);

[sphereXi,sphereYi,sphereZi]=ind2sub(size(voxelRadii),find(voxelRadii<rad_mm));      % matrix indices of voxels within the sphere

relRoi=[sphereXi-ceil(size(xrange,2)/2), sphereYi-ceil(size(yrange,2)/2), sphereZi-ceil(size(zrange,2)/2)];    % RELATIVE matrix indices of voxels within the sphere

end%function
