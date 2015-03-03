function RDMs=squareRDMs(RDMs_ltv)
% converts set of row-vector RDMs_ltv to square form (despite being
% rows, RDMs are stacked along the 3rd dimension, just as square RDMs
% would be. this avoids ambiguity when the RDMs_ltv is square and could
% be either a single RDM or a number of vectorized RDMs.)
% RDMs may be bare or wrapped with meta-data in a struct object. they
% will be returned in the same format as passed.
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

if isstruct(RDMs_ltv)
    % wrapped
    RDMs_ltv_struct=RDMs_ltv;
    RDMs_ltv=unwrapRDMs(RDMs_ltv_struct);
    
    nRDMs=size(RDMs_ltv,3);
    RDMs=[];
    for RDMI=1:nRDMs
        RDMs=cat(3,RDMs,squareRDM(RDMs_ltv(:,:,RDMI)));
    end
    
    RDMs=wrapRDMs(RDMs,RDMs_ltv_struct);
else
    % bare
    nRDMs=size(RDMs_ltv,3);
    RDMs=[];
    for RDMI=1:nRDMs
        RDMs=cat(3,RDMs,squareform(vectorizeRDM(RDMs_ltv(:,:,RDMI))));
    end
end

end%function
