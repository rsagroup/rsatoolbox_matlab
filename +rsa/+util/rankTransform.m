function rankArray=rankTransform(array,scale01)

% transforms the array 'array' by replacing each element by its rank in the
% distribution of all its elements. setting scale01 would scale the
% elements to [0-1] after ranktrasforming all the non-NaN entries.
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

if ~exist('scale01','var'), scale01=false; end;

nonNan_LOG=~isnan(array);
set=array(nonNan_LOG); % column vector

[sortedSet, sortedIs]=sort(set);

rankArray=nan(size(array));
nonNan_IND=find(nonNan_LOG);
rankArray(nonNan_IND(sortedIs))=1:numel(set);

if scale01
    rankArray=rankArray/numel(set);
end

end%function
