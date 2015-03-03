function RGBim=mat2RGBimage(mat,cols,clims)

% Given a matrix mat, a colormap cols, and colorscale-limiting values
% clims, this function returns an RGB image RGBim in which the mat is
% displayed using the colormap, such that the first color in cols maps onto
% matrix value clims(1) and the last color in cols maps onto clims(2), with
% a linear scaling in between these values. Matrix values smaller than
% clim(1) are represented by cols(1,:) and matrix values larger than
% clims(2) by cols(end,:).

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

nCols=size(cols,1);
linRange=linspace(clims(1),clims(2),nCols);
RGBim=interp1(linRange,cols,mat);
smallerLOG=repmat(mat<=clims(1),[1 1 3]);
RGBim(smallerLOG)=[cols(1,1)*ones(sum(smallerLOG(:))/3,1);cols(1,2)*ones(sum(smallerLOG(:))/3,1);cols(1,3)*ones(sum(smallerLOG(:))/3,1)];
largerLOG=repmat(mat>=clims(2),[1 1 3]);
RGBim(largerLOG) =[cols(end,1)*ones(sum(largerLOG(:))/3,1);cols(end,2)*ones(sum(largerLOG(:))/3,1);cols(end,3)*ones(sum(largerLOG(:))/3,1)];

end%function
