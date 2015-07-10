function imageRDMs(RDMs,varargin)
% imageRDMs(RDMs,varargin);
%
% Visualizes one or many RDMs, and uses the *.RDM fields for each structure
% which could be sqaure or vectorized
%
% 
% Joern Diedrichsen 2015
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*;
import rsa.fig.*;
import rsa.util.*;

% define default behavior
Opt.figureNumber   = gcf;       % Overprint the current figure
Opt.clims          = [];        % Color scale
Opt.transformFcn   = 'ssqrt';   % Transformfunction
Opt.showColorbar   = true;
Opt.aspect         = 2/3;
Opt.imagelabels    = [];
Opt.colourScheme   = []; 
Opt = rsa.getUserOptions(varargin,Opt);

% Deal with different input arguments

allMin = inf;
allMax = -inf;

nRDMs   = length(RDMs); 
RDMs    = rsa.rdm.vectorizeRDMs(RDMs);
allRDMs = vertcat(RDMs.RDM);
switch (Opt.transformFcn)
    case ''
        % Nothing to do
    case 'rank'
        scale01(rankTransform_equalsStayEqual(allRDMs,2))
    otherwise
        allRDMs=feval(Opt.transformFcn,allRDMs);
end;


if isempty(Opt.clims)
    Opt.clims = [min(allRDMs(:)) max(allRDMs(:))];
end;

% display dissimilarity matrices
h=figure(Opt.figureNumber);
set(h,'Color','w');

[nVerPan nHorPan]=paneling(nRDMs+1,Opt.aspect);
clf;

for i=1:nRDMs
    subplot(nVerPan,nHorPan,i);
    thisRDM = rsa.rdm.squareRDM(allRDMs(i,:));
    
    % Determine alpha data to make nans invisible
    alpha = ~isnan(thisRDM);
    
    image(thisRDM,'CDataMapping','scaled','AlphaData',alpha);
    set(gca,'CLim',Opt.clims,'CLimMode','manual');
    if (~isempty(Opt.colourScheme))
        colormap(gca, Opt.colourScheme);
    end; 
    
    set(gca,'XTick',[],'YTick',[]);
    
    if isstruct(RDMs) && isfield(RDMs,'name')
        title(['\bf' deunderscore(RDMs(i).name)]);
    else 
        title(['\bf' sprintf('%d',i')]);
    end;
    axis square off;
    
    % If image labels, add them
    if ~isempty(Opt.imagelabels)
        addImageSequenceToAxes(gca,Opt.imagelabels);
    end
end


% add color bar
if Opt.showColorbar
    subplot(nVerPan,nHorPan,nRDMs+1); cla;
    imagesc(thisRDM,Opt.clims);  
    cla;
    set(gca,'XLim',[0 10],'YLim',[0 10]); 
    ht=text(5,5,{['\bfdissimilarity matrices'],'not rank-transformed'},'HorizontalAlignment','Center','FontUnits','normalized');
    set(ht,'FontSize',.06);
    axis square off;
    colorbar;
end
