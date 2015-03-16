function h=pageFigure(figI,paperSizeORheightToWidth,proportionOfScreenArea,horPos0123,landscapeFig)
% creates a page figure with figure with specific properties.
% h=pageFigure([figI,paperSizeORheightToWidth,proportionOfScreenArea,horPos1234]);
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

if ~exist('figI'), newFigure=gcf; figI=newFigure.Number; end
if ~exist('paperSizeORheightToWidth')||isempty(paperSizeORheightToWidth), paperSizeORheightToWidth='A4'; end
if ~exist('proportionOfScreenArea')||isempty(proportionOfScreenArea), proportionOfScreenArea=0.25; end
if ~exist('landscapeFig')||isempty(landscapeFig), landscapeFig=false; end

if figI
    h=figure(figI(1));
else
    h=figure;
end

if ~exist('horPos0123')||isempty(horPos0123), horPos0123=mod(mod(h.Number,10),4); end

set(h,'Color','w');
%set(h,'WindowStyle','docked');

if isstr(paperSizeORheightToWidth)
    if strcmp(paperSizeORheightToWidth,'A4')
        heightToWidth=sqrt(2)/1;
    elseif strcmp(paperSizeORheightToWidth,'legal')
        heightToWidth=14/8.5;
    elseif strcmp(paperSizeORheightToWidth,'letter')
        heightToWidth=11/8.5;
    end
else
    heightToWidth=paperSizeORheightToWidth;
end

if landscapeFig
    heightToWidth=1/heightToWidth;
end

lbwh = get(0,'ScreenSize');
screenArea=lbwh(3)*lbwh(4);

figWidth=sqrt(screenArea*proportionOfScreenArea/heightToWidth);
figHeight=heightToWidth*figWidth;

left=figWidth*horPos0123;
% left=(lbwh(3)-figWidth)/2;
bottom=(lbwh(4)-figHeight)/2;

set(h,'Position',[left bottom figWidth figHeight])
% [left, bottom, width, height]

set(h,'PaperPositionMode','auto'); % 'auto' here prevents resizing when the figure is printed.

end%function
