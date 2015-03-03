function RDMs_squareNbare=stripNsquareRDMs(RDMs_clad)
% strips and squares a set of RDMs embedded in a structure also specifying
% names etc. to return the bare RDMs stacked along the 3rd (which are
% assumed to be in field RDM) or, if RDMs are already bare, they are
% returned as passed. 
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

if isstruct(RDMs_clad)
    RDMs_squareNbare=[];
    nRDMs=numel(RDMs_clad);
    for RDMI=1:nRDMs
        RDMs_squareNbare=cat(3,RDMs_squareNbare,squareRDM(RDMs_clad(RDMI).RDM));
    end
else
    RDMs_squareNbare=squareRDMs(RDMs_clad);
end