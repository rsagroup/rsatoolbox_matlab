function addHeading(heading,figI,x,y)
% adds a heading to a figure. The user also has the option of specifying
% the location of the inserted heading text.

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

% replace underscores
if iscell(heading)
    for lineI=1:numel(heading)
        line=heading{lineI};
        line(line==95)='-';
        heading{lineI}=line;
    end
else
    heading(heading==95)='-';
end    

if ~exist('figI','var')
    newFigure = gcf;
    figI=newFigure.Number;
end
pageFigure(figI);

if ~exist('x','var'), x=1.11; end
if ~exist('y','var'), y=1.08; end

h=axes('Parent',gcf); hold on;
set(h,'Visible','off');
axis([0 1 0 1]);

% add heading(s)
text(x,y,heading,'HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12,'FontWeight','bold','Color','k');

end%function
