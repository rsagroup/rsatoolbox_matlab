function roi=mask2roi(mask)
% gives the oordinates of the marked points in an input matrix  mask.
% this function assumes that any non-zero entries are part of the roi.
% For example when the input is a 3D binary mask, the output would be
% triplets containing the coordinates of the points (in matrix space) where
% the binary mask is set to one.
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

[x,y,z]=ind2sub(size(mask),find(mask));
roi=[x,y,z];

end%function
