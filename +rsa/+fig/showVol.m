function showVol(volORmap, title, figI, right, skipNslices, nHorPanels, nVerPanels, labels, labelPos)

% USAGE:        showVol(vol, [title, figI, right, skipNslices, nHorPanels, nVerPanels, labels, labelPos])
%
% FUNCTION:     to display a volume as a set of slices
%
% ARGUMENTS:
% vol           the volume as a stack of true-color slices:
%               X by Y by 3 by Z. the third dimension encodes 
%               the color component (red, green, blue).
%               anatomically, the X axis points to the left,
%               the Y axis to the back of the brain, and the
%               Z axis up.
%
% [title]       string to be used as a title (optional)
%
% [figI]        figure number (optional), generates autonumbered next
%               figure if this argumant is missing or empty.
%
% [right]       a string specifying the orientation of the
%               slices. if right is 'right', then right is 
%               right, i.e. the view onto the slices is from
%               above. if right is 'left' or simply 'wrong',
%               then right is left, i.e. the view onto the 
%               slices is from below.
%
% [skipNslices] optional number of slices to be skipped between slices to
%               be shown. defaults to 0 (all slices shown).
%
% labels        structured array of string labels (optional)
%
% labelPos      nLabels by 3 matrix of label positions (optional)

% DEBUG
% vol(1,1,:,1)=[1 0 0];
% vol(10,1,:,1)=[1 0 0];
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

if ndims(volORmap)==3 % if its a map...
    vol=map2vol(volORmap);
else % must be a vol
    vol=volORmap;
end

if ~exist('skipNslices','var')
    sliceIs=1:size(vol,4);
else
    sliceIs=1:skipNslices+1:size(vol,4);
end

vol=vol(:,:,:,sliceIs);

if nargin<7 || nHorPanels==0
    nHorPanels=ceil(sqrt(size(vol,4)));
    nVerPanels=ceil(sqrt(size(vol,4)));
end

if ~exist('right','var')
    right='right';
end

if ~exist('figI','var') || (exist('figI','var') && isempty(figI))
    figI=0;
end

if ~exist('title','var')
    title='';
end

if figI(1)
    figure(figI(1));
    if numel(figI)>1
        subplot(figI(2),figI(3),figI(4)); cla;
        if numel(figI)>4
            % vol=vol(:,:,:,figI(5));
            vol=vol(5:75,40:75,:,figI(5));
        end
    else
        clf;
    end
else
    figI=figure;
end

set(figI(1),'Color','w');
%set(figI,'WindowStyle','docked');

margin=0;
width=(1-(nHorPanels-1)*margin)/nHorPanels;
height=(1-(nVerPanels-1)*margin)/nVerPanels;

[sx,sy,sc,sz]=size(vol);

for sliceI=1:size(vol,4)
    %subplot(nVerPanels,nHorPanels,sliceI);
    left=mod(sliceI-1,nHorPanels)*(width+margin);
    bottom=1-ceil(sliceI/nHorPanels)*(height+margin);
    
    if numel(figI)<2
        subplot('Position',[left bottom width height]);
    end
    
    if strcmp(right,'right')
        im=flipdim(permute(vol(:,:,:,sliceI),[2 1 3]),2);
    else
        im=permute(vol(:,:,:,sliceI),[2 1 3]);
    end
    image(im);
    
    axis equal;
    axis([0 size(im,2)-1 0 size(im,1)-1])

    set(gca,'Visible','off');
end

if exist('labels','var') && strcmp(right,'right')
    labelPos(:,1)=sx-labelPos(:,1); % flip x
end

h=axes('Parent',gcf); hold on;
set(h,'Visible','off');
axis([0 1 0 1]);

text(0.5,-0.04,['Right is ',right,'.'],'HorizontalAlignment','Center','Color',[0.7 0.7 0.7],'FontSize',7);
% if strcmp(right,'right')
%     text(0.5,-0.04,'Right is right.','HorizontalAlignment','Center','Color',[0.7 0.7 0.7],'FontSize',7);
% else
%     text(0.5,-0.04,'Right is left.','HorizontalAlignment','Center','Color',[0.7 0.7 0.7],'FontSize',7);
% end


% text(0.5,-0.08,title,'HorizontalAlignment','Center','FontSize',7,'FontWeight','bold','Color',[.4 .4 .4]);
text(0.5,1,title,'HorizontalAlignment','Center','FontSize',14,'FontWeight','bold','Color',[.8 .8 .8]);
%text(0.5,1.05,title,'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','Color','k');

if exist('labels','var')
    for labelI=1:length(labels);
        sliceI=labelPos(labelI,3);
        left=mod(sliceI-1,nHorPanels)*(width+margin);
        bottom=1-ceil(sliceI/nHorPanels)*(height+margin);
        subplot('Position',[left bottom width height]);
        
        % text(labelPos(labelI,1),labelPos(labelI,2),labels(labelI),'BackgroundColor','w','HorizontalAlignment','center');
        text(labelPos(labelI,1),labelPos(labelI,2),labels(labelI),'Color','w','FontWeight','bold','HorizontalAlignment','center');
        text(labelPos(labelI,1),labelPos(labelI,2),labels(labelI),'Color','k','FontWeight','normal','HorizontalAlignment','center');
    end
end

end%function
