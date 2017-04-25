function vol=addBinaryMapToVol(vol, map, col)

% USAGE:        vol=addBinaryMapToVol(vol, map, col)
%
% FUNCTION:     to superimpose the binary map map to the true-color
%               volume vol. where map==1 the corresponding voxel of volume
%               vol is marked in the color col. all other voxels retain
%               their color.
%
% ARGUMENTS:
% vol           the volume as a stack of true-color slices: X by Y by 3 by Z.
%               the third dimension encodes the color component (red, green,
%               blue).
%               anatomically, the X axis points to the left, the Y axis to
%               the back of the brain, and the Z axis up.
%
% map           a statistical map as a 3D array of double-precision floats
%
% col           triple row vector [R G B] specifying the color, in which voxels
%               that contain a 1 in the map map are marked in the volume vol.
%               values should be element of [0,1].
%
% RETURN VALUE: true-color volume (X by Y by 3 by Z) with the binary
%               map superimposed in color col.
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

%% MARK THE VOXELS
toBeMarked_INDs=find(map==1);

volXYZ3=permute(vol,[1 2 4 3]);

  redmap=volXYZ3(:,:,:,1);
greenmap=volXYZ3(:,:,:,2);
 bluemap=volXYZ3(:,:,:,3);

  redmap(toBeMarked_INDs)=col(1);
greenmap(toBeMarked_INDs)=col(2);
 bluemap(toBeMarked_INDs)=col(3);

volXYZ3(:,:,:,1)=  redmap;
volXYZ3(:,:,:,2)=greenmap;
volXYZ3(:,:,:,3)= bluemap;

vol=permute(volXYZ3,[1 2 4 3]);

end%function
