function plotDotsWithTextLabels(coords_xy, localOptions)

% plots the set text labels in cell array 'textLabels' at the positions
% 'coords_xy' in the current plot in font size 'fontSize' in black (by
% default) or using category colors specified by optional further arguments
% 'categoryColumns', 'categoryColors', and 'categoryLabels'.
%
% localOptions
%	textLabels: cell of strings
%	dotColours: nItems x 3 matrix
%	fontSize: int
%	categoryColumns
%	categoryColours
%	categoryLabels
%
% CW 5-2010: Added dots ability
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

localOptions = setIfUnset(localOptions, 'dotSize', 8);

nItems=size(coords_xy,1);

if ~exist('localOptions', 'var') || isempty(localOptions.textLabels)
	for i = 1:nItems
		localOptions.textLabels{i} = num2str(i);
	end%for:i
end%if

if ~isfield(localOptions, 'fontSize')
        localOptions.fontSize=11;
end%if:fontSize

hold on;

for i = 1:numel(localOptions.textLabels)
	localOptions.textLabels{i} = deunderscore(localOptions.textLabels{i});
end%for:i

if isfield(localOptions, 'categoryColumns')

    nCategories=size(localOptions.categoryColumns,2);
    if ~isfield(localOptions, 'categoryColors')
        localOptions.categoryColors=randomColor(nCategories);
	end%if:categoryColours

    for categoryI=1:nCategories

        cCategoryItemIs=find(logical(localOptions.categoryColumns(:,categoryI)));

        for itemI=cCategoryItemIs
            text(coords_xy(itemI,1),coords_xy(itemI,2),...
                localOptions.textLabels{itemI},'Color',localOptions.categoryColors(categoryI,:),'FontSize',localOptions.fontSize);
		end%for:itemI
	end%for:categoryI

    if isfield(localOptions, 'categoryLabels')
        xlim=get(gca,'XLim');
        ylim=get(gca,'YLim');
        yr=max(range(ylim),range(xlim)/4);
        for categoryI=1:nCategories
           text(xlim(1)+(categoryI-1+.5)*range(xlim)/nCategories,ylim(1)-0.1*yr,localOptions.categoryLabels(categoryI),'FontWeight','bold','Color',localOptions.categoryColors(categoryI,:),'HorizontalAlignment','Center');
		end%for:categoryI
	end%if:categoryLabels
else
    % plot text labels in black
    for itemI=1:nItems
        text(coords_xy(itemI,1),coords_xy(itemI,2),...
            localOptions.textLabels{itemI},'Color','k','FontSize',localOptions.fontSize);
	end%for:itemI
end%if:categoryColumns

if isfield(localOptions, 'dotColours')
	for itemI = 1:nItems
		plot(coords_xy(itemI,1), coords_xy(itemI,2),'o','MarkerFaceColor',localOptions.dotColours(itemI, :),'MarkerEdgeColor','none','MarkerSize', localOptions.dotSize);
	end%for:itemI
end%if:dotColours
axis off

end%function
