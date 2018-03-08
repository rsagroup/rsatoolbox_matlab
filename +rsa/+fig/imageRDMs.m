function imageRDMs(RDMs,varargin)
% imageRDMs(RDMs,varargin);
%
% Visualizes one or many RDMs, and uses the *.RDM fields for each structure
% can deal with a number of different input formats:
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
Opt.singleRDM      = 0;         % If singleRDM flag is set, it only plots the current axis 
Opt.lines          = []; 
Opt = rsa.getUserOptions(varargin,Opt);

allMin = inf;
allMax = -inf;

% Deal with different input formats 
if (isstruct(RDMs)) 
    nRDMs   = length(RDMs); 
    if (nRDMs>1)            % this is the struct array with one struct per RDM 
        RDMs = rsa.util.struct2dataframe(RDMs); 
    elseif (nRDMs==1) 
        nRDMs = size(RDMs.RDM,1); 
    end; 
else 
    nRDMs = size(RDMs,1); 
    a=RDMs; 
    RDMs=[]; 
    RDMs.RDM=a;
end; 

switch (Opt.transformFcn)
    case ''
        % Nothing to do
    case 'rank'
        scale01(rankTransform_equalsStayEqual(RDMs.RDM,2))
    otherwise
        RDMs.RDM=feval(Opt.transformFcn,RDMs.RDM);
end;


if isempty(Opt.clims)
    Opt.clims = [min(RDMs.RDM(:)) max(RDMs.RDM(:))];
end;

if (nRDMs>1 || Opt.singleRDM==0)
    % Make new Figure     
    h=figure(Opt.figureNumber);
    set(h,'Color','w');
    [nVerPan nHorPan]=paneling(nRDMs+1,Opt.aspect);
    clf;
end; 

for i=1:nRDMs
    if (nRDMs>1 || Opt.singleRDM==0)
        subplot(nVerPan,nHorPan,i);
    end; 
    thisRDM = rsa.rdm.squareRDM(RDMs.RDM(i,:));
    
    % Determine alpha data to make nans invisible
    alpha = ~isnan(thisRDM);
    
    image(thisRDM,'CDataMapping','scaled','AlphaData',alpha);
    set(gca,'CLim',Opt.clims,'CLimMode','manual');
    if (~isempty(Opt.colourScheme))
        colormap(gca, Opt.colourScheme);
    end; 
    
    set(gca,'XTick',[],'YTick',[]);
    
    % Draw lines 
    if (~isempty(Opt.lines)) 
        drawline(Opt.lines); 
        drawline(Opt.lines,'dir','horz'); 
    end; 
    
    if isfield(RDMs,'name')
        if (iscell(RDMs.name))
            title(['\bf' deunderscore(RDMs.name{i})]);
        else 
            title(['\bf' deunderscore(RDMs.name)]);
        end; 
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
if Opt.showColorbar && Opt.singleRDM==0
    subplot(nVerPan,nHorPan,nRDMs+1); cla;
    imagesc(thisRDM,Opt.clims);  
    cla;
    set(gca,'XLim',[0 10],'YLim',[0 10]); 
    ht=text(5,5,{['\bfdissimilarity matrices'],'not rank-transformed'},'HorizontalAlignment','Center','FontUnits','normalized');
    set(ht,'FontSize',.06);
    axis square off;
    colorbar;
end