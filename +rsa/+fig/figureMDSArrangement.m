function figureMDSArrangement(RDM, userOptions, localOptions)

% FUNCTION
%       draws a multidimensional scaling (MDS) arrangement reflecting the
%       dissimilarity structure of the items whose
%       representational dissimilarity matrix is passed in argument RDM.
%       the function can draw the MDS result as an arrangement of text
%       labels (default), colored dots, or icons (e.g. the experimental
%       stimuli or, more generally, icons denoting the experimental
%       conditions). which ones of these visualizations are produced is
%       controlled by setting the fields of the options argument.
%
%		Written by Niko Kriegeskorte, edited by Cai Wingfield
%
% USAGE
%       showMDSarrangement(RDM, userOptions, localOptions)
%
% ARGUMENTS
% RDM
%       the items' dissimilarity matrix as a struct RDM,
%       or in square-matrix or upper-triangular-vector format.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.criterion
%                        The criterion which will be minimised to optimise the
%                        MDS arrangement. Defaults to metric stress.
%
% localOptions
%       optional struct whose fields control which visualizations are
%       produced and provide the required additional information. the
%       fields (like the struct as a whole) are optional. if options is not
%       passed or none of the fields are set, the MDS arrangement will be
%       drawn using text labels, where each label is a number indexing an
%       activity pattern in the order defined by the RDM.
%
%       the the key fields are:
%       [figI_textLabels] figure index for visualization as text-label
%                       arrangement. defaults to 1.
%       [figI_catCols]  figure index for visualization as colored-dot
%                       arrangement, where the colors code for category.
%                       this visualization is omitted if this argument is
%                       missing.
%       [figI_icons]    figure index for visualization as an icon
%                       arrangement. this visualization is omitted if this 
%                       argument is missing.
%       [figI_shepardPlots] figure index for shepard plot (a scatterplot
%                       relating the dissimilarities in the original space
%                       to the distances in the 2S MDS arrangement. 
%       please note: if the figI_... arguments are quadruples instead of
%                       single figure indices, then the last three numbers
%                       specify the subplot. see selectPlot.m for details.
%       [MDScriterion]  the cost function minimized by the MDS. the default
%                       value is 'metricstress'. see documentation of
%                       mdscale.m for details.
%       [rubberbandGraphPlot] boolean argument. if this is set to 'true',
%                       then a rubberband graph is plotted in the
%                       background to visualize the distortions incurred by
%                       the dimensionality reduction (for details, see
%                       kriegeskorte et al. 2008).
%
%       further fields needed for visualization as category-color-coded
%       arrangement of either dots or text labels:
%       [contrasts]     a matrix, whose columns define categories of 
%                       items as index vectors (column height
%                       == number of items, 1 indicates
%                       present, 0 indicates absent). (these category
%                       definitions are a special case of a general-linear-
%                       model contrast -- hence the name of this argument.)
%       [categoryIs]    list of intergers referring to columns of argument 
%                       'contrasts' and thereby selecting which categories 
%                       of items are to be included in the 
%                       visualization.
%       [categoryColors]nCategories-by-3 matrix, whose rows define the
%                       colors as RGB triples. there is one row per
%                       category.
%       [categoryLabels] text labels for the categories. if these are
%                       provided a legend will show what the color-coded
%                       categories are.
%
%       further field needed for visualization as icon arrangement:
%       [icons]         the icon images for visualization as an icon
%                       arrangement.

%
% EDITS
%         3-2010: CW: It now works with userOptions and localOptions, and some
%                     of the options fieldnames have been changed.
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

%% define defaults
appendFlag = 0;
% warning('off', 'IgnoringExtraEntries');

if ~exist('localOptions','var')||isempty(localOptions), localOptions=struct; end
localOptions=setIfUnset(localOptions,'figI_textLabels',1);   

if isstruct(RDM)
    RDMname=RDM.name;
else
    RDMname='[unnamed dissimilarity matrix]';
end
description{1}=[RDMname ', ' userOptions.criterion];

figIs=[];

%% perform multidimensional scaling (MDS)
%D = unwrapSimmats(squareSimmats(RDM));
D = unwrapRDMs(squareRDMs(RDM));

nDims=2;

try
    [pats_mds_2D, stress, disparities] = mdscale(D, nDims,'criterion',userOptions.criterion, 'options', struct('MaxIter', 1000));
catch
    try
        [pats_mds_2D, stress, disparities] = mdscale(D, nDims,'criterion','stress');
		% pats_mds_2D is an [nPoints 2]-sized matrix of coordinates
        description{1}=[RDMname ', reverted to stress: ' userOptions.criterion ' failed)'];
    catch
        try
            D2=D+0.2;
            D2(logical(eye(length(D)))) = 0;
            [pats_mds_2D, stress, disparities] = mdscale(D2, nDims,'criterion',userOptions.criterion);
            description{1}=[RDMname ', ' userOptions.criterion ' , added 0.2 to distances to avoid colocalization'];
        catch
            description{1}=[RDMname ', MDS failed...'];
        end
    end   
end    


if isfield(localOptions,'figI_shepardPlots')&&~isempty(localOptions.figI_shepardPlots)
    figIs=[figIs,localOptions.figI_shepardPlots];
    shepardPlot(D,disparities,pdist(pats_mds_2D),localOptions.figI_shepardPlots,['\fontsize{12}MDS (' description{1} ')']);
end


%% plot MDS arrangement using text labels
if isfield(localOptions,'figI_textLabels')&&~isempty(localOptions.figI_textLabels)
    
    figIs=[figIs,localOptions.figI_textLabels(1)]; [hf,ha]=selectPlot(localOptions.figI_textLabels);

    % plot rubberband plot in the background
    if ~isfield(localOptions,'rubberbandGraphPlot')
        if size(D,2)<10
            localOptions.rubberbandGraphPlot=1;
        else
            localOptions.rubberbandGraphPlot=0;
        end
    end

    if localOptions.rubberbandGraphPlot
        rubberbandGraphPlot(pats_mds_2D,D);
    end

    if ~exist('contrasts','var')
        % categories undefined: plot all text labels in black
		veryLocalOptions.textLabels = localOptions.textLabels;
		veryLocalOptions.dotColours = localOptions.dotColours;
		if isfield(userOptions, 'dotSize')
			veryLocalOptions.dotSize = userOptions.dotSize;
		end%if
        plotDotsWithTextLabels(pats_mds_2D,veryLocalOptions);
		
		% Plot convex hulls
		if isfield(localOptions, 'convexHulls')
		
			hold on;
		
			categories = unique(localOptions.convexHulls);
			nCategories = numel(categories);
			for hullI = 1:nCategories
				
				thisHullLabel = categories(hullI);
				withinHullPointIs = find(localOptions.convexHulls == thisHullLabel);
				withinHullPointCoords = pats_mds_2D(withinHullPointIs, :);
				hullBoundaryPointIs = convhull(withinHullPointCoords(:,1), withinHullPointCoords(:,2));
				plot(withinHullPointCoords(hullBoundaryPointIs, 1), withinHullPointCoords(hullBoundaryPointIs, 2), 'LineWidth', 2, 'Color', localOptions.dotColours(withinHullPointIs(1), :));
				
			end%for:hullI
			
			hold off;
			
		end%if
		
        axis equal off;
    else
        % categories defined: plot text labels in category colors
        nCategories=size(localOptions.contrasts,2);
		jetCols=jet;
        localOptions=setIfUnset(localOptions,'categoryColors',jetCols(round(linspace(1,64,nCategories)),:));

        defaultLabels=cell(nCategories,1);
        for i=1:nCategories
            defaultLabels(i)={num2str(i)};
        end
        localOptions=setIfUnset(localOptions,'categoryLabels',defaultLabels);

        figIs=[figIs,pageFigure(localOptions.figI_catCols(1))]; [hf,ha]=selectPlot(localOptions.figI_catCols);

        if ~isfield(localOptions,'categoryIs')
            localOptions.categoryIs=1:size(localOptions.contrasts,2);
        end

	if isfield(localOptions, 'fontSize')
		fontSize = localOptions.fontSize;
	else
		fontSize = 11;
	end%if:fontSize

        plotTextLabels(pats_mds_2D,localOptions.textLabels,fontSize,localOptions.contrasts(:,localOptions.categoryIs),localOptions.categoryColors(localOptions.categoryIs,:),localOptions.categoryLabels(localOptions.categoryIs));
        axis equal off;
        if isstruct(RDM), title(['\bf',deunderscore(RDM.name)]); end
end

axis([min(pats_mds_2D(:,1)) max(pats_mds_2D(:,1)) min(pats_mds_2D(:,2)) max(pats_mds_2D(:,2))]);

if isfield(localOptions, 'titleString')
	title(['\fontsize{12}' localOptions.titleString]);
end%if

% Handle the figure
if isfield(userOptions, 'colourmap')
		set(gcf,'Colormap', userOptions.colourmap);
	end%if

	% Then export and/or close figures appropriately
	gotoDir(userOptions.rootPath);
	fileName = [userOptions.analysisName '_' localOptions.fileName];
	handleCurrentFigure(fileName, userOptions);
	clear thisFileName
end

end%function


% %% plot MDS arrangement using category-color-coded dots
% if isfield(options,'figI_catCols')&&~isempty(options.figI_catCols) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nCategories=size(options.contrasts,2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     options=setIfUnset(options,'categoryColors',jetCols(round(linspace(1,64,nCategories)),:)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     defaultLabels=cell(nCategories,1);
%     for i=1:nCategories
%         defaultLabels(i)={num2str(i)};
%     end
%     options=setIfUnset(options,'categoryLabels',defaultLabels); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     figIs=[figIs,pageFigure(options.figI_catCols(1))]; [hf,ha]=selectPlot(options.figI_catCols);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     if ~isfield(options,'categoryIs')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         options.categoryIs=1:size(options.contrasts,2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
%     
%     if exist('pats_mds_2D','var')
%         plotDots(pats_mds_2D,options.contrasts(:,options.categoryIs),options.categoryColors(options.categoryIs,:),options.categoryLabels(options.categoryIs)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         axis tight equal off;
%         if isstruct(RDM), title(['\bf',deunderscore(RDM.name)]); end
%     end
% end
% 
% 
% %% plot MDS arrangement using icons
% if isfield(options,'figI_icons')&&~isempty(options.figI_icons) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     figIs=[figIs,pageFigure(options.figI_icons(1))]; [hf,ha]=selectPlot(options.figI_icons); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     if ~isfield(userOptions,'rubberbandGraphPlot')
%         if size(D,2)<10
%             userOptions.rubberbandGraphPlot=1;
%         else
%             userOptions.rubberbandGraphPlot=0;
%         end
%     end
%     
%     if userOptions.rubberbandGraphPlot
%         rubberbandGraphPlot(pats_mds_2D,D);
%     end
% 
%     drawImageArrangement(options.icons,pats_mds_2D); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     axis tight equal off;
% 
% end


% warning('on', 'IgnoringExtraEntries');
