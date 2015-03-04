function drawImageArrangement(imageStruct,coords_xy,imageAreaProportion,transparentCol)

% FUNCTION
%       to draw an arrangement of images (imageStruct) into the current
%       axes at specified locations (coords_xy).
%
% USAGE
%       drawImageArrangement(imageStruct,coords_xy[,imageAreaProportion,transparentCol])
%
% ARGUMENTS
% imageStruct
%       structured array containing the images in the field 'image' (in
%       matlab's 'CData' format for the built-in function 'image', e.g. 3D
%       arrays, horizontal-by-vertical-by-RGB). for the moment, the images
%       are assumed to be square.
%
% coords_xy
%       matrix of image coordinates. each row contains the x and y
%       coordinates for one image with the row indices corresponding to the
%       indices of imageStruct. 
%
% [imageAreaProportion]
%       optional argument that defaults to 1. if this is a singleton, it
%       controls the area of the images relative to the area of the entire
%       arrangement (i.e. a rectangular bounding box around coords_xy). 1
%       indicates that the sum of the areas of all images equals the area
%       of the arrangement. if this argument has two elements, then the
%       first element is ignored and the second one directly specifies the
%       image width.
%
% [transparentCol]
%       optional argument that defaults to [128 128 128 2]. the first three
%       elements specify a color (RGB) to be rendered as transparent. if
%       this argument has a fourth element, then the transparent color will
%       fade smoothly from opaque to transparent around the part of the
%       image that has other colors, according to a gaussian whose sigma is
%       transparentCol(4).
%
% STATE OF DEVELOPMENT
%       draft stage
%       (should be generalized further to allow better control of drop shadows,
%       arbitrary alpha channels, and nonsquare images)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% preparations
if ~exist('imageAreaProportion','var'), imageAreaProportion=1; end
if ~exist('transparentCol','var'), transparentCol=[128 128 128 2]; end

if numel(imageStruct)~=size(coords_xy,1)
    warning('drawImageArrangement:  numel(imageStruct)~=size(coords_xy,1)');
    %pause;
end


nImages=size(coords_xy,1);

axis equal;


%% compute image size
% images assumed to be square
if numel(imageAreaProportion)<2
    if nImages>1
        boundingBoxArea=max(prod(range(coords_xy)),max(range(coords_xy))^2/10);
    else
        xlim=get(gca,'XLim');
        ylim=get(gca,'YLim');
        boundingBoxArea=range(xlim)*range(ylim);
    end
    totalImageArea=boundingBoxArea*imageAreaProportion;
    imageWidth=sqrt(totalImageArea/nImages);
else
    imageWidth=imageAreaProportion(2);
end

%% arrange the images
if numel(transparentCol)==4
    % smooth alpha channel
    hsize=5*transparentCol(4);
    sigma=1*transparentCol(4);
    kernel=fspecial('gaussian', hsize, sigma);
end

for imageI=1:nImages

    %[xs,ys,rgb3]=size(imageStruct(imageI).image);

    transparent=imageStruct(imageI).image(:,:,1)==transparentCol(1) & imageStruct(imageI).image(:,:,2)==transparentCol(2) & imageStruct(imageI).image(:,:,3)==transparentCol(3);

    if numel(transparentCol)==4
        % smooth alpha channel
        opacity=imfilter(double(1-transparent),kernel);
    else
        opacity=~transparent;
    end   
    
    image('CData',imageStruct(imageI).image,'XData',[coords_xy(imageI,1)-imageWidth/2, coords_xy(imageI,1)+imageWidth/2],'YData',[coords_xy(imageI,2)+imageWidth/2, coords_xy(imageI,2)-imageWidth/2],'AlphaData',opacity);
end

% axis tight equal off;
% set(gca,'xtick',[],'ytick',[]);
