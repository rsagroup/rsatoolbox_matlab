function F = setPapertoFigPos(F)
% This function change's the paper size of a figure to match
% the space the figure takes up on the page. This is mainly
% since it allows you to set the 'position' argument to the
% desired size and then have the resulting ps/eps/pdfs come out
% without a bunch of whitespace.
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

if nargin < 1
	F = gcf;
end

%paperunits and units tend to be different (cm and pixels)
set(F,'units',get(F,'paperunits'));

% we want to scale the papersize to match the figure's size
figpos = get(F,'paperposition');
set(F,'papersize',figpos(3:4));

% Of course now the figure's paperposition is probably outside
% the paper. We know the paper is now the same size as the fig
% so by setting the distances to 0 we should now be smack in the
% middle (?)
% If matlab's 'position' setting works, it should sort itself out
set(F,'paperposition',[0 0 figpos(3:4)]);

end%function
