function plotTextLabels(coords_xy,textLabels,fontSize,categoryColumns,categoryColors,categoryLabels)

% plots the set text labels in cell array 'textLabels' at the positions
% 'coords_xy' in the current plot in font size 'fontSize' in black (by
% default) or using category colors specified by optional further arguments
% 'categoryColumns', 'categoryColors', and 'categoryLabels'.
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

nItems=size(coords_xy,1);

if ~exist('textLabels','var')||isempty(textLabels)
    for itemI=1:nItems
        textLabels{itemI}=num2str(itemI);
    end
end



if ~exist('fontSize','var')||isempty(fontSize)
        fontSize=15;
end

hold on;

for i = 1:numel(textLabels)
	textLabels{i} = deunderscore(textLabels{i});
end

if exist('categoryColumns','var')

    nCategories=size(categoryColumns,2);
    if ~exist('categoryColors','var')||isempty(categoryColors)
        categoryColors=randomColor(nCategories);
    end

    for categoryI=1:nCategories

        cCategoryItemIs=find(logical(categoryColumns(:,categoryI)));

        for itemI=cCategoryItemIs
            text(coords_xy(itemI,1),coords_xy(itemI,2),...
                textLabels{itemI},'Color',categoryColors(categoryI,:),'FontSize',fontSize);
        end
    end

    if exist('categoryLabels','var')
        xlim=get(gca,'XLim');
        ylim=get(gca,'YLim');
        yr=max(range(ylim),range(xlim)/4);
        for categoryI=1:nCategories
            % vertical
            % text(xlim(1)+0.05*range(xlim),ylim(1)+categoryI*range(ylim)/(nCategories+1),categoryLabels(categoryI),'FontWeight','bold','Color',categoryColors(categoryI,:));

            % horizontal
            text(xlim(1)+(categoryI-1+.5)*range(xlim)/nCategories,ylim(1)-0.1*yr,categoryLabels(categoryI),'FontWeight','bold','Color',categoryColors(categoryI,:),'HorizontalAlignment','Center');
        end
        % legend(categoryLabels,'Location','NorthWest');
    end
else % ~exist('categoryColumns','var')
    % plot text albels in black
    for itemI=1:nItems
        text(coords_xy(itemI,1),coords_xy(itemI,2),...
            textLabels{itemI},'Color','k','FontSize',fontSize);
    end
end

end%function
