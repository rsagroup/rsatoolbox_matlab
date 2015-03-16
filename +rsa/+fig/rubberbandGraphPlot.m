function rubberbandGraphPlot(coords_xy,distmat,connectionAreaProportionORcolorCoding)
% given a set of coordinates (coords_xy, two columns one for x and another
% for y coordinates) and the corresponding distance matrix (distmat) and a
% ratio for scaling the areas (connectionAreaProportionORcolorCoding,
% greater values corresponding to larger areas), this function creates a
% rubberband graph in which points defined by cords_xy are connected by
% rubberbands for which the area is proportinal to the exact dissimilarity.
% This feature enables displaying the exact dissimilarities in an mds plot
% that may be distorted due to dimensionality reduction (see also
% shepardPlot.m)
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

%% control variables
nDistortionBins=11;

%% preparations
if exist('connectionAreaProportionORcolorCoding','var')
    if strcmp(class(connectionAreaProportionORcolorCoding),'double')
        connectionAreaProportion=connectionAreaProportionORcolorCoding;
        colorCoding=false;
    elseif strcmp(class(connectionAreaProportionORcolorCoding),'char') && strcmp(connectionAreaProportionORcolorCoding,'color')
        colorCoding=true;
        connectionAreaProportion=.5;
    end
else
    connectionAreaProportion=.5;
    colorCoding=false;
end
distmat_utv=vectorizeRDM(distmat);
nPoints=size(coords_xy,1);

% pixelsPerInch=get(0,'ScreenPixelsPerInch');
% pointsPerInch=72; % definition of 'point' unit
% pixelsPerPoint=pixelsPerInch/pointsPerInch;

hold on;
plot(coords_xy(:,1),coords_xy(:,2),'LineStyle','none');
axis tight equal;
set(gca,'Units','points'); lbwh_pts=get(gca,'Position');
unitsPerPoint=range(get(gca,'XLim'))/lbwh_pts(3);


%% compute 2D distmat
distmat2D_utv=pdist(coords_xy,'euclidean');


%% compute connection thicknesses
boundingBoxArea=prod(range(coords_xy));
areas_utv=distmat_utv/sum(distmat_utv)*boundingBoxArea*connectionAreaProportion;
thickness_utv=areas_utv./distmat2D_utv; % area=distmat2D*thickness, proportional to distmat


%% map from upper-triangular-vector (utv) form of distance matrix to point indices i, j
[i,j]=ndgrid(1:nPoints,1:nPoints);
i_utv=vectorizeRDM(i);
j_utv=vectorizeRDM(j);


%% draw the connections
if colorCoding
    squeezedCol=[0 .5 1];
    stretchedCol=[1 0 0];
    undistortedCol=[.9 .9 .9];
    distortionCols=interp1(1:3,[squeezedCol;undistortedCol;stretchedCol],linspace(1,3,nDistortionBins));
    lw=boundingBoxArea/sum(distmat2D_utv)/unitsPerPoint*connectionAreaProportion;
else
    lc=[.9 .9 .9]; % line color
end

mnt=min(thickness_utv);
mxt=max(thickness_utv);
thicknessCtrs=linspace(mnt,mxt,nDistortionBins);
binWidth=thicknessCtrs(2)-thicknessCtrs(1);
thicknessEdges=[thicknessCtrs-binWidth/2,thicknessCtrs(end)+binWidth/2];

if mxt-mnt<1e-6 % if no distortion...
    % select coords of all pairs
    sourcePositions_xy=coords_xy(i_utv,:);
    targetPositions_xy=coords_xy(j_utv,:);
    Z=[repmat(-1,size(sourcePositions_xy(:,1)')); repmat(-1,size(sourcePositions_xy(:,1)'))];

    % draw the lines
    if colorCoding
        line([sourcePositions_xy(:,1)'; targetPositions_xy(:,1)'],[sourcePositions_xy(:,2)'; targetPositions_xy(:,2)'],Z,'Color',undistortedCol,'LineWidth',lw);
    else
        thickness_pts=mnt/unitsPerPoint;
        line([sourcePositions_xy(:,1)'; targetPositions_xy(:,1)'],[sourcePositions_xy(:,2)'; targetPositions_xy(:,2)'],Z,'Color',lc,'LineWidth',thickness_pts);
    end
else % if distortions present...
    for distortionBinI=1:nDistortionBins

        % find pairs to be connected by lines of this thickness
        pairs_utvLOG=thicknessEdges(distortionBinI)<thickness_utv & thickness_utv<=thicknessEdges(distortionBinI+1);
        thickness_pts=mean(thicknessEdges(distortionBinI:distortionBinI+1))/unitsPerPoint;

        is=i_utv(pairs_utvLOG);
        js=j_utv(pairs_utvLOG);

        % select coords of those pairs
        sourcePositions_xy=coords_xy(is,:);
        targetPositions_xy=coords_xy(js,:);
        Z=[repmat(-1,size(sourcePositions_xy(:,1)')); repmat(-1,size(sourcePositions_xy(:,1)'))];

        % draw the lines
        if colorCoding
            line([sourcePositions_xy(:,1)'; targetPositions_xy(:,1)'],[sourcePositions_xy(:,2)'; targetPositions_xy(:,2)'],Z,'Color',distortionCols(distortionBinI,:),'LineWidth',lw);
        else
            line([sourcePositions_xy(:,1)'; targetPositions_xy(:,1)'],[sourcePositions_xy(:,2)'; targetPositions_xy(:,2)'],Z,'Color',lc,'LineWidth',thickness_pts);
        end
    end
end


% sourceTargetVecs=targetPositions_xy-sourcePositions_xy;
%     % draw lines as patches
%     % draw the line
%     corners=[sourcePos+(null(lineVec)*(lineWidth/2))';
%         sourcePos-(null(lineVec)*(lineWidth/2))';
%         targetPos-(null(lineVec)*(lineWidth/2))';
%         targetPos+(null(lineVec)*(lineWidth/2))'];
% 
%     patch(corners(:,1),corners(:,2),lineCol,'EdgeColor','none'); % function patch draws one polygon per column

end%function
