function showVoxObj(mapORcoords, figI, subplotTriple, adjustMargin, axisLabels, materialDef)

% USAGE
%       showVoxObj(mapORcoords[, figI, subplotTriple, adjustMargin, axisLabels])
%
% FUNCTION
%       to display a binary voxel object defined by mapORcoords
% 
% ARGUMENTS
% mapORcoords 
%       this can be a 3D binary map or a list of 3D coordinates
%       (nVoxels by 3).
%
% [adjustMargin]
%       should be left unspecified or set to 1 (default) for normal
%       display. the margin is then adjusted to the bounding box of the
%       object (expanded or reduced), such that there is a 1 voxel margin
%       around the object before upsampling.
%
% [axisLabels]
%       optional structure with fields 'x','y' and 'z', which have string
%       values defining the axis labels.
%       

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

objCol=[0 1 0];
backgroundCol=[1 1 1];

if nargin==1
    figI=1248; clf;
    subplotTriple=[1 1 1];
end

if nargin<4
    adjustMargin=1;
end

if size(mapORcoords,2)==3 && ndims(mapORcoords)==2
    % convert the list of coordinates to a map
    map=zeros(max(mapORcoords(:,1)),max(mapORcoords(:,2)),max(mapORcoords(:,3)));
    indices=sub2ind(size(map),mapORcoords(:,1),mapORcoords(:,2),mapORcoords(:,3));
    map(indices)=1;
else
    map=mapORcoords;
    clear mapORcoords;
end

if ~exist('materialDef','var'), materialDef=[0   0.5   2.3  1 0.4]; end

if adjustMargin
    [x,y,z]=ind2sub(size(map),find(map));

    minX=min(x); maxX=max(x);
    minY=min(y); maxY=max(y);
    minZ=min(z); maxZ=max(z);

    sizeX=maxX-minX+1;
    sizeY=maxY-minY+1;
    sizeZ=maxZ-minZ+1;

    mapWithMargin=zeros(sizeX+2,sizeY+2,sizeZ+2);
    mapWithMargin(2:1+sizeX,2:1+sizeY,2:1+sizeZ)=map(minX:maxX,minY:maxY,minZ:maxZ);
else
    mapWithMargin=map;
end
    
h=figure(figI);
set(h,'Color',backgroundCol);

subplot(subplotTriple(1),subplotTriple(2),subplotTriple(3));
cla;

upsampleFactor=4;
largerMapWithMargin=upSample(mapWithMargin,upsampleFactor);
p = patch(isosurface(permute(largerMapWithMargin,[2 1 3]),0.5));
%isonormals(x,y,z,v,p)
set(p,'FaceColor',objCol,'EdgeColor','none');

axis equal;
%axis off;
axmargin=0;
axis([-axmargin size(largerMapWithMargin,1)+axmargin -axmargin size(largerMapWithMargin,2)+axmargin -axmargin size(largerMapWithMargin,3)+axmargin]);
camva(10);
if exist('axisLabels','var')&&~isempty(axisLabels)
    xlabel(axisLabels.y);
    ylabel(axisLabels.x);
    zlabel(axisLabels.z);
else
    xlabel('y');
    ylabel('x');
    zlabel('z');    
end

%light('Position',[1 0 0],'Style','infinite');
%set(h,'FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5)
view(160,30);

% view(94,4);
% %lightangle(5,4)

camlookat(p);
camlight('headlight');
camlight('left');
material(materialDef);
lighting phong;
title([num2str(length(find(map))),' voxels']);

disp(['This object has ',num2str(length(find(map))),' voxels.']);

end%function
