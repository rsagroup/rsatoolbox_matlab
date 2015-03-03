function rankMat=rankTransform_equalsStayEqual(mat,scale01)

% transforms the matrix mat by replacing each element by its rank in the
% distribution of all its elements. if scale01 is set the ranked elements
% would be scaled between 0 to 1. One property of this function, as
% indicated by its name, is that equal values would stay equal after
% ranking and scaling. This feature is not present in the other function of
% the toolbox (rankTransform.m). The way the function does this is by
% assigning the mean scaled rank to all the equal entries of the input
% matrix (mat). NaNs are ignored in this process.
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

if ~exist('scale01','var'), scale01=1; end;

nonNan_LOG=~isnan(mat);
set=mat(nonNan_LOG);

[sortedSet, sortedIs]=sort(set);

rankMat=nan(size(mat));
nonNan_IND=find(nonNan_LOG);
rankMat(nonNan_IND(sortedIs))=1:numel(set);

if scale01==1
    % scale into [0,1]
    rankMat=(rankMat-1)/(numel(set)-1);
elseif scale01==2
    % scale into ]0,1[
    % (best representation of a uniform distribution between 0 and 1)
    rankMat=(rankMat-.5)/numel(set);
end


%% fix artefactual differences
uniqueValues=unique(mat(nonNan_LOG));

for uniqueValueI=1:numel(uniqueValues)
    cValueEntries_LOG=mat==uniqueValues(uniqueValueI);
    rankMat(cValueEntries_LOG)=mean(rankMat(cValueEntries_LOG));
end

end%function
