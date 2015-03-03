function figureDendrogram(RDM, userOptions, localOptions)

% Generates a dendrogram representation of an RDM and handles the figure
% according to user preferences.
%
% figureDendrogram(RDM, userOptions, localOptions)
%
%        RDM --- The distance matrix basis for the dendrogram.
%                RDM should be in squareform: a [nConditions nConditions]-sized
%                matrix.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.conditionColours
%                        A [nConditions 3]-sized matrix indicating an [R G B]
%                        tripple colour for each condition.
%                userOptions.saveFigurePDF
%                        A boolean value. If true, the figure is saved as a PDF.
%                        Defaults to false.
%                userOptions.saveFigurePS
%                        A boolean value. If true, the figure is saved as a PS.
%                        Defaults to false.
%                userOptions.saveFigureFig
%                        A boolean value. If true, the figure is saved as a
%                        MATLAB .fig file. Defaults to false.
%                userOptions.displayFigures
%                        A boolean value. If true, the figure remains open after
%                        it is created. Defaults to true.
%
%        localOptions --- Further options.
%                localOptions.linkageType
%                        Defaults to 'single'.
%                localOptions.colorThreshold
%                        Defaults to uncoloured.
%
% Figures may be exported according to preferences.
%
% Cai Wingfield 3-2010
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

if ~isfield(localOptions, 'linkageType')
	localOptions.linkageType = 'single';
end

[hf ha figI] = selectPlot(localOptions.figureNumber);
if isfield(localOptions, 'colorThreshold')
	[H_ignore T_ignore labelReordering] = dendrogram(linkage(vectorizeRDM(RDM), localOptions.linkageType), 0, 'labels', localOptions.labels, 'colorThreshold', localOptions.colorThreshold, 'orientation', 'left');
else
	[H_ignore T_ignore labelReordering] = dendrogram(linkage(vectorizeRDM(RDM), localOptions.linkageType), 0, 'labels', localOptions.labels, 'orientation', 'left');
end%if

if isfield(userOptions, 'conditionColours');
	hold on;
	if size(userOptions.conditionColours, 1) == size(squareRDM(RDM), 1)
		x = xlim(gca);
		for condition = 1:size(squareRDM(RDM), 1)
			plot(x(1), labelReordering(condition), 'o', 'MarkerFaceColor', userOptions.conditionColours(condition, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
		end%for:condition
	else
		error('figureDendrogram:wrongNumberOfConditionColours', 'Number of conditions is not the same as the number colours of userOptions.conditionColours. Please fix it!');
	end%if
end%if

if isfield(localOptions, 'titleString')
	title(['\bf' localOptions.titleString]);
end%if

% Then export and/or close figures appropriately
gotoDir(userOptions.rootPath);% , 'Figures'
fileName = [userOptions.analysisName '_' localOptions.fileName];
handleCurrentFigure(fileName, userOptions);
clear thisFileName

end%function
