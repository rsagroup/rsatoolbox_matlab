function [hf,ha,figI]=selectPlot(figI)
% - activates or creates figure figI(1)
% - activates or creates plot [figI(2:4)],
%   or [1 1 1] if figI(2:4) are missing
% - sets background color to white
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

if ~exist('figI','var'); figI=0; end

if figI
    hf=figure(figI(1));
    if numel(figI)==4;
        ha=subplot(figI(2),figI(3),figI(4));
    else
        ha=subplot(1,1,1);
    end
else
    hf=figure;
    figI(1)=hf;
    ha=subplot(1,1,1);
end
set(hf,'Color','w');  

end%function
