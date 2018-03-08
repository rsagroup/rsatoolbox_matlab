function g=gaussian_nk(FWHM, nSamples, x, mu, mass)
% USAGE
%       g=gaussian_nk(FWHM[, nSamples, x, mu, mass])
%
% FUNCTION
%       returns a row-vector gaussian of full width at half maximum FWHM.
%       nSamples is chosen automatically unless it is passed as an
%       argument. if the other params are unspecified the domain is 
%                       -(nSampels-1)/2:(nSampels-1)/2 
%       and nSamples must be odd, so there is a central sample for the
%       maximum of the gaussian. the center mu is then placed at zero and
%       the mass==sum(g) is set to 1. the unit of the FWHM is assumed to be
%       the sampling period (i.e. the interval between adjacent samples),
%       unless the domain x is passed. if the domain x is passed, it needs
%       to have nSamples many points and it will define the units of FWHM. 
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

if ~exist('nSamples','var')
    nSamples=FWHM*4+1;
    % this automatic choice makes sense when the FWHM because the gaussian
    % will be very low this far out (<0.0001). the unit of the FWHM is assumed to be
    % the sampling period, i.e. the interval between adjacent samples.
end    
    
nLatSamples=(nSamples-1)/2;
if ~exist('x','var')
    if nLatSamples~=floor(nLatSamples) % nSamples even -> no central sample for the maximum
        error('ERROR in gaussian: nSamples even -> no central sample for the maximum');
    end
    x=-nLatSamples:nLatSamples;
end

if ~exist('mu','var')
    mu=0;
end

if ~exist('mass','var')
    mass=1;
end



s=1/2.3548; % relation between standard deviation (in case of pdf) and fwhm

g=exp(-x.^2/(2*(FWHM*s)^2));
g=g/sum(g)*mass;

% figure(200); clf; plot(x,g,'r','LineWidth',3)

end%function
