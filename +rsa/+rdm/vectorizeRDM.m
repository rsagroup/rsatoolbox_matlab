function RDM=vectorizeRDM(RDM)
% converts RDM to lower-triangular form (row vector) (or leaves it in that
% form) in cases where the input(RDM) is a vector, the output would be the
% same as the input (row or column vector). if the input is a 2D matrix
% (e.g. an RDM) the output would be a (row) vector of the lower-triangular
% entries.
% this function does not accept a wrapped RDM as its input.
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

if size(RDM,1)==size(RDM,2)
    RDM(logical(eye(size(RDM))))=0; % fix diagonal: zero by definition
    RDM=squareform(RDM);
end

end%function
