function  [binRDM, nCatCrossingsRDM]=categoricalRDM(categoryVectors,figI,monitor)

% Given a single vector of category indices, this function returns a binary
% representational dissimilarity matrix (binRDM) containing a zero for each
% pair of conditions falling into the same category, and a one for each
% pair of conditions falling into different categories.
%
% Given a set of column vectors that define categories, this function
% returns a binary RDM indicating the condition pairs straddling at least
% one category boundary. A 0 indicates that the dissimilarity at that
% location is between two conditions that are within the same category on
% all category vectors. A 1 indicates that the dissimilarity is between two
% conditions that span separate categories according to at least one of the
% category vectors.
% 
% The additional return value nCatCrossingsRDM contains a count of the
% number of category boundaries that divide each pair of conditions. This
% number ranges from 0 (both conditions in the same category according to
% all category vectors) to m, the width of category vectors (the conditions
% are in different categories according to all category vectors).
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
if ~exist('monitor','var'), monitor=true; end
if ~exist('figI','var'), figI=500; end
if min(size(categoryVectors))==1, categoryVectors=categoryVectors(:); end
[nCond nCats]=size(categoryVectors);


%% count category crossings for each pair of conditions
nCatCrossingsRDM=zeros(nCond,nCond);

for catI=1:nCats
    cCatBinRDM=repmat(categoryVectors(:,catI),[1 nCond])~=repmat(categoryVectors(:,catI)',[nCond 1]);
    nCatCrossingsRDM=nCatCrossingsRDM+cCatBinRDM;
end

binRDM=logical(nCatCrossingsRDM);


%% visualise
if monitor
    showRDMs(concatRDMs(binRDM,nCatCrossingsRDM),figI);
end


end%function
