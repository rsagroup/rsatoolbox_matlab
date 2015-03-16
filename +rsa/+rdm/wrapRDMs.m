function RDMs_struct=wrapRDMs(RDMs,refRDMs_struct)
% wraps similarity matrices RDMs (in square or upper triangle form)
% into a structured array with meta data copied from refRDMs_struct
% (which needs to have the same number of RDMs).(if they are already
% wrapped then the wrapping (metadata) is replaced by that of
% refRDMs_struct.
% generally in cases where one wants to wrap a number of dissimilarity
% vectors or matrices, they can define them as different 'RDM' fields of an
% input structure. In cases where the user wants to define names or specify
% colours for the dissimilarity vector or matrices they neeed to specify
% this in the refRDMs_struct, otherwise generic names and values would be
% selected and returned in the 'name' and 'color' fields.
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

if isstruct(RDMs)
    % wrapped already, but replace the wrapping
    nRDMs=numel(RDMs);
else
    nRDMs=size(RDMs,3);
end    

if ~exist('refRDMs_struct','var')
    for RDMI=1:nRDMs
        refRDMs_struct(RDMI).name='[unnamed RDM]';
        refRDMs_struct(RDMI).color=[0 0 0];
    end
end

RDMs_struct=refRDMs_struct;
if isstruct(RDMs)
    % wrapped already, but replace the wrapping
    for RDMI=1:nRDMs
        RDMs_struct(RDMI).RDM=RDMs(RDMI).RDM;
    end
else
    % RDMs need wrapping
    for RDMI=1:nRDMs
        RDMs_struct(RDMI).RDM=RDMs(:,:,RDMI);
    end
end

end%function
