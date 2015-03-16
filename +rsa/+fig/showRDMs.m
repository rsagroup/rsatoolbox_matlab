function showRDMs(RDMs,figI,rankTransform01,clims,showColorbar, aspect, imagelabels, colourScheme)
% visualizes one or many RDMs. RDMs is a structure of RDMs (wrapped RDMs).
% The 'RDM' and 'name' subfileds are required for this function only. figI
% is the figure number in which the RDMs are displayed. if rankTransform01
% is defined, the dissimilarities would be rankTransformed before
% visualization. clims is a 1x2 vector specfying the lower and upper limits
% for displayed dissimilarities. if showColorbar is set, the colorbar would
% also be displayed in a separate panel. aspect: number of horizontal
% panels/number of vertical panels. colourScheme: the colour scheme for
% displaying dissimilarities (defined by RDMcolormap by default). 
% imagelabels: a structure containg the RGB values for images that would be
% displayed on an RDM. 
% CW: 5-2010, 7-2010. HN: 9-2013
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% define default behavior
if ~exist('figI','var'), figI=500; end
if ~exist('clims','var'), clims=[]; end
if ~exist('rankTransform01','var'), rankTransform01=true; clims=[0 1]; end
if ~exist('showColorbar','var'), showColorbar=true; end
if ~exist('aspect', 'var') || isempty(aspect), aspect = 2/3; end

%% handle RDM types

colourScheme = RDMcolormap;
if isstruct(RDMs)
	[rawRDMs nRDMs] = unwrapRDMs(RDMs);
else
	rawRDMs = RDMs;
	nRDMs = size(rawRDMs, 3);
end%if:isstruct(RDMs)

% rawRDMs is now [nC nC nR] or [1 t(nC) nR]

cellRDMs = cell(1, nRDMs);

allMin = inf;
allMax = -inf;

for RDMi = 1:nRDMs

	thisRDM = rawRDMs(:, :, RDMi);

	if max(size(thisRDM)) == numel(thisRDM)
		% Then it's in ltv form
		% So we've got a symmetric RDM
		% Square it
		thisRDM = squareform(thisRDM); % squareform leaves 0s on the diagonal, which is what we want
		RDMtype{RDMi} = '';
	else
		% Then it's a square RDM
		if isequalwithequalnans(thisRDM, thisRDM')
			% It's symmetric
			if ~any(diag(thisRDM)) && ~any(isnan(diag(thisRDM))) % 0s on the diagonal
				% It's a regular RDM
				RDMtype{RDMi} = '';
			elseif isnan(thisRDM(1,1)) % (1,1) is a nan
				% It's a very lucky cv2RDM
				RDMtype{RDMi} = '[cv2RDM] ';
			else
				% symmetric, not all 0s or all nans on diagonal
				warning(['RDM ' num2str(RDMi) ' is symmetric, but has neither all 0s or all NaNs on the diagonal.  Not sure how to deal with this.']);
                RDMtype{RDMi} = '[anomalous diag.] '
			end
		else
			% It's not symmetric
			if isnan(thisRDM(1,1)) % (1,1) is a nan
				% It's a cv2RDM
				RDMtype{RDMi} = '[cv2RDM] ';
			else
				% It's a sdRDM
				RDMtype{RDMi} = '[sdRDM] ';
			end
		end
	end

	% Work out extreme values
	offDiagonalEntries = thisRDM(eye(size(thisRDM, 1))==0);
	thisMin = nanmin(offDiagonalEntries);
	thisMax = nanmax(offDiagonalEntries);

	allMin = min(allMin, thisMin);
	allMax = max(allMax, thisMax);

	cellRDMs{1,RDMi} = thisRDM;

end%for:RDMi

% cellRDMs is now a [1 nRDMs]-sized cell of square RDMs (of various kinds and sizes)

if isempty(clims)
	clims = [allMin allMax];
end%if:clims

%% display dissimilarity matrices
h=figure(figI(1)); set(h,'Color','w');

if numel(figI)<4
    [nVerPan nHorPan]=paneling(nRDMs+1,aspect);
    subplotOffset=0;
    clf;
else
    nVerPan=figI(2);
    nHorPan=figI(3);
    subplotOffset=figI(4)-1;
end
    
for RDMi=1:nRDMs

	subplot(nVerPan,nHorPan,RDMi+subplotOffset); cla;

	thisRDM = cellRDMs{RDMi};

	% Determine alpha data to make nans invisible
	alpha = ~isnan(thisRDM);

	if rankTransform01
		image(scale01(rankTransform_equalsStayEqual(thisRDM,1)),'CDataMapping','scaled','AlphaData',alpha);
		set(gca,'CLim',[0 1],'CLimMode','manual');
	else
		image(thisRDM,'CDataMapping','scaled','AlphaData',alpha);
		set(gca,'CLim',clims,'CLimMode','manual');
	end
	
	if exist('colourScheme', 'var')
		colormap(gca, colourScheme);
	end%if
	
	set(gca,'XTick',[],'YTick',[]);
	
	if isstruct(RDMs)
		title(['\bf' RDMtype{RDMi} deunderscore(RDMs(RDMi).name)]);
	elseif ~isempty(RDMtype{RDMi})
		title(['\bf' RDMtype{RDMi}]);
	end;
	axis square off;
	
	% If image labels, add them
	if exist('imagelabels','var') && ~isempty(imagelabels)
		addImageSequenceToAxes(gca,imagelabels);
	end
	
end


%% add color bar
if showColorbar
    subplot(nVerPan,nHorPan,nRDMs+1+subplotOffset); cla;
	tempRDM = cellRDMs{1,1};
	n = size(tempRDM,1);
    if rankTransform01
        imagesc(rankTransform(squareRDM(tempRDM),1),clims);  cla;
        ht=text(n/2,n/2,{['\bfeach dissimilarity matrix (',num2str(n),'^2)'], 'separately rank-transformed', 'and scaled into [0,1]'},'HorizontalAlignment','Center','FontUnits','normalized');
	else
        imagesc(tempRDM,clims);  cla;
        ht=text(n/2,n/2,{['\bfdissimilarity matrices (',num2str(n),'^2)'],'not rank-transformed'},'HorizontalAlignment','Center','FontUnits','normalized');
    end
    set(ht,'FontSize',.06);
    axis square off;
    %colormapJet4Print;
    colorbar;
end

end%function
