function Y=temporallySmoothTimeSpaceMatrix(Y,FWHM)
% Y=temporallySmoothTimeSpaceMatrix(Y,FWHM)
% temporally smooths the time-by-space matrix Y with a gaussian kernel of
% full width at half maximum FWHM
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

%% define smoothing kernel
kernel=gaussian_nk(FWHM)';


%% smooth temporally
Y=conv2(Y,kernel,'same');

end%function
