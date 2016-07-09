function r=raeSpearmanCorr(x,y,toleratedStandardError)

% This function computes the "random-among-equals" Spearman correlation
% (RAE-Spearman r), which is the Pearson correlation applied to ranks with
% random ranks assigned to sets of entries in either variable that have
% equal values. In case there are no ties, the RAE-Spearman correlation is
% equal to the conventional Spearman correlation. If there are ties, then
% the RAE-Spearman correlation, uses a random ranks within each equality
% class (rather than average ranks as in the Spearman correlation).
%
% This modification is motivated by the scenario of comparing model
% predictions to data, where some models predict ties. The Spearman (and
% also the Pearson and Kendall tau b and tau c) grants an advantage to
% models predicting ties over models predicting at chance or slightly
% better than chance for the same pairs of entries. As a result, these
% correlation coefficients grant a higher score to models that avoid
% detailed predictions by predicting ties over models that make more
% accurate and precise predictions. A simplified model can beat the true
% model when both are compared to noisy data. The Kendall tau a does not
% have this drawback, but it has time complexity O(n^2), where n is the
% number of entries -- making it ill-suited for extensive permutation or
% bootstrap testing of complex model predictions.
%
% The RAE-Spearman, like the Kendall tau a, penalises the prediction of
% ties by replacing the ties with random guesses. Like the tau a it rewards
% a model (in terms of the expected value of the r) for predicting the
% direction of a pair of entries whenever it has *any information* about
% the pair of entries (i.e. whenever the probability of a correct guess is
% larger than 0.5). Both the RAE-Spearman correlation and Kendall's tau a
% penalise ties as much as guessing.
%
% The RAE-Spearman has the advantage of being faster to compute than the
% Kendall tau a. Its time complexity is O(n*log(n)) (because sorting is
% O(n*log(n)) and the Pearson correlation is O(n)). It has the disadvantage
% of being nondeterministic, because random ranks are chosen. The parameter
% nReps offers a simple remedy, controlling how many times the RAE-Spearman
% is computed and results averaged. 
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

if ~exist('toleratedStandardError','var'), toleratedStandardError=0.001; end

x=x(:); y=y(:);
validEntryIs = ~isnan(x)&~isnan(y);
x=x(validEntryIs); y=y(validEntryIs);

nEstimates_max=10;
rs=nan(nEstimates_max,1);
nEstimatesPerRound=3;
nEstimates=0;
cStandardError=inf;

while nEstimates<nEstimates_max && toleratedStandardError<cStandardError
    
    for estimateI=1:nEstimatesPerRound
        xy_rank=tiedrankrandom([x y]);
        rs(nEstimates+estimateI)=fastPearsonCorr(xy_rank(:,1),xy_rank(:,2));
    end
    nEstimates=nEstimates+nEstimatesPerRound;
    
    cStandardError=std(rs(1:nEstimates))/sqrt(nEstimates);
end

r=mean(rs(1:nEstimates));

end%function
