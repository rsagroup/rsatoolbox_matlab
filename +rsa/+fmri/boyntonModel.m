function hirf=boyntonModel(res)

%returns a vector containing the hemodynamic impulse response function
%as estimated by boynton, engel, glover and heeger (1996) for two subjects' V1.
%the parameter res controls the temporal resolution (time bin width is 1 ms).

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

delta=2500; %as in BV (averages of the subjects?)
tau=1250;
n=3;


t=0:res:20000;

hirf=(t/tau).^(n-1).*exp(-t/tau)/(tau*factorial(n-1));
hirf=[zeros(1,floor(delta/res)),hirf];

end%function
