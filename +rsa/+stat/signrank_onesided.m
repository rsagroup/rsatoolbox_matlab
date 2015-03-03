function [p,h] = signrank_onesided(x);
% returns the p-value from a one-sided Wilcoxon signed rank test comparing
% against zero. The exact method is chosen to get the p-values with greater
% resolution.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

[p_signrank,h] = signrank(x,[],'alpha',0.05,'method','exact');
if median(x) > 0
     p_signrank = p_signrank/2;
else
    p_signrank = 1-p_signrank/2;
end

p = p_signrank;

end%function
