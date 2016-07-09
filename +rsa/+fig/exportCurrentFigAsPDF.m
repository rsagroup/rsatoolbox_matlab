function exportCurrentFigAsPDF(filespec,userOptions)
% exports the current figures to the file [filespec,'.pdf'] in pdf
% format. 
% filespec is the full path (including the file name) of the file to be
% saved.
% userOptions can be used to modify the dots per inch (dpi) or the
% tightness of the figure in the printed page.
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

if ~isfield(userOptions,'dpi'), userOptions.dpi = 300; end;
if ~isfield(userOptions,'tightInset'), userOptions.tightInset = false; end;

% If we are trimming the paper size to match the figure
if userOptions.tightInset
    setPapertoFigPos;
end
appendFlag = 0;
switch appendFlag
    case 0
        print('-dpdf',sprintf('-r%d',userOptions.dpi),filespec);
    case 1
        print('-dpdf',sprintf('-r%d',userOptions.dpi),'-append',filespec);
    case 3
        print('-dpdf',sprintf('-r%d',userOptions.dpi),[filespec,'_',num2str(gcf)]);
        print('-dpdf',sprintf('-r%d',userOptions.dpi),'-append',filespec);
end

end%function
