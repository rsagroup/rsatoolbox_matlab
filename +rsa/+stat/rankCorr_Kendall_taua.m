function taua=rankCorr_Kendall_taua(a,b)
% computes the Kendall's tau a correlation coefficient between the input
% vectors (a and b). NaN entries would be removed from both.
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

%% preparations
a=a(:);b=b(:);
validEntryIs = ~isnan(a)&~isnan(b);
a=a(validEntryIs);b=b(validEntryIs);
n=size(a,1);


%% compute Kendall rank correlation coefficient tau-a
K = 0;
for k = 1:n-1
    pairRelations_a=sign(a(k)-a(k+1:n));
    pairRelations_b=sign(b(k)-b(k+1:n));
    K = K + sum(pairRelations_a.*pairRelations_b);
end
taua=K/(n*(n-1)/2); % normalise by the total number of pairs 

end%function
