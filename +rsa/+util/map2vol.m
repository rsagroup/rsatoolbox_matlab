function vol=map2vol(map)
%
% converts a map (e.g. statistical t-map) to a true-color volume. 
% This function can be helpful for displaying a map using the other
% function called showVol.m
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

if strcmp(class(map),'single') || strcmp(class(map),'double')
    % scale into RGB range [0,1] for display
    map=map-min(map(:)); % move min to zero
    map=map/max(map(:)); % scale max to one
    vol=double(permute(repmat(map,[1 1 1 3]),[1 2 4 3]));
else
    % make vol same class as map (e.g. logical)
    vol=permute(repmat(map,[1 1 1 3]),[1 2 4 3]);
end

end%function
