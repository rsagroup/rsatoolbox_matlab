function largerMap=upSample(map,factor)

% USAGE:        largerMap=upSample(map,factor)
%
% FUNCTION:     to upsample the 3D matrix map by repeating each value
%               factor times
%
% AUTHOR:       nk
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

mapDim=size(map);
largerMapDim=factor*mapDim;
if numel(largerMapDim)==2
    largerMapDim=[largerMapDim,1];
end

[lxI, lyI, lzI]=meshgrid(1:largerMapDim(1),1:largerMapDim(2),1:largerMapDim(3));
lxI=ceil(lxI/factor);
lyI=ceil(lyI/factor);
lzI=ceil(lzI/factor);

indices=sub2ind(mapDim,lxI, lyI, lzI);

largerMap=map(indices);

end%function
