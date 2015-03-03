function [] = image_thr(im,thr);
% this function takes an image (im) and a threshold (thr) as its input. It
% thresholds the image and displays the values *bellow* the threshold in
% red. The other values are shown in grayscale.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

nCols = 256;clims = [0 1];
cols = colorScale([0 0 0;1 1 1],nCols);
RGBim=mat2RGBimage(im,cols,clims);
[u,v] = find(im <= thr);
for i=1:numel(u)
    RGBim(u(i),v(i),:) = [1 0 0];
end
image(RGBim);
colormap(cols);axis square

end%function
