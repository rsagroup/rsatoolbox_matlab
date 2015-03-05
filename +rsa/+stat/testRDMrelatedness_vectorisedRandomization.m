%	testRDMrelatedness_vectorisedRandomization
%
%	USAGE
%
%	p = ...
%	[r p] = ...
%	[r p p_conv] = ...
%
%		= testRDMrelatedness_vectorisedRandomization(RDMa,RDMb[, correlationType, nPermutations])
%
%	EXPLANATION
%
%	This function tests the null hypothesis that similarity matrices A and
%	B are unrelated. the test simulates a null distribution of correlations
%	between A and B by means of randomization of the conditions labels
%	of the similarity matrix.
%
%	Completely based on Russell Thompson's genius method of vectorising the
%	permutation test.  (Really.  I just changed the variable names and
%	beefed up the comments.)
%
%	Cai Wingfield 4-2010

function [varargout]=testRDMrelatedness_vectorisedRandomization(RDMa,RDMb,varargin)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% Deal with varargin

if numel(varargin) == 0
	correlationType = 'Spearman';
	nPermutations = 10000;
elseif numel(varargin) == 1
	correlationType = varargin{1};
	nPermutations = 10000;
elseif numel(varargin) == 2
	correlationType = varargin{1};
	nPermutations = varargin{2};
end%if:nargin

%% Initial tests and variables setup

if size(RDMa,1) ~= size(RDMa,2)
	RDMa = squareform(RDMa);
end%if
if size(RDMb,1) ~= size(RDMb,2)
	RDMb = squareform(RDMb);
end%if

nConditions = size(RDMa, 1);

if size(RDMb, 1) ~= nConditions
	error('Both input RDMs must be of the same size.');
end%if

%% %%%%%%%%%%%%%%%%%%%%%% %%
%% The randomisation test %%
%% %%%%%%%%%%%%%%%%%%%%%% %%

%% Create the permuted labels

% Create a [nPermutations nConditions]-matrix where each row is a new random permutation
[ignore, permutationLables] = sort(rand(nPermutations,nConditions),2);

% Add the identity permutation on top, then sort the order of the remaining permutations lexicographically
permutationLables = unique([1:nConditions;permutationLables],'rows');

% Convert these labels into indices referencing elements of RDMa:
% The repmat takes the permutationsLables and copies them nCondition times next to eachother
% The kron makes a (similarly sized) matrix containing values which, when added to the previous matrix, convert the row indices to cell indices approprite for permuted rows AND colums (clever!!)
referenceIndices = repmat(permutationLables, 1, nConditions) + kron((permutationLables - 1) * nConditions, ones(1, nConditions));

% Now get the unique cells
uniqueCells = tril(reshape(1:nConditions^2,nConditions,nConditions),-1);
uniqueCells = uniqueCells(uniqueCells~=0);
	
% Pick out the lower-triangular cells from RDMb
RDMb = RDMb(uniqueCells);

% Pick out ALL the permuted, lower-triangular cells from RDMa
RDMaRand = RDMa(referenceIndices(:,uniqueCells)');
RDMa = RDMa(uniqueCells);

% Correlate RDMa and RDMb
[rows, columns]=size(RDMaRand);
switch correlationType
	case {'Correlation', 'Pearson'}
		RDMaRand = sum((RDMaRand-repmat(mean(RDMaRand),rows,1)).*repmat(RDMb-mean(RDMb),1,columns))./((std(RDMaRand).*std(RDMb)).*(rows-1));
	case {'Spearman', 'Rank'}
		[RDMaRand(:), ignore] = tiedrank(RDMaRand(:));
		RDMaRand = sum((RDMaRand-repmat(mean(RDMaRand),rows,1)).*repmat(RDMb-mean(RDMb),1,columns))./((std(RDMaRand).*std(RDMb)).*(rows-1));
end%switch:correlationType

% Position of real correlation in permuted distribution
[ignore, referenceIndices] = sort(RDMaRand);
p = max([1-(find(referenceIndices == 1) / columns) 1 / columns]);

%% %%%%%%%% %%
%% The rest %%
%% %%%%%%%% %%

%% Deal with varargout

if nargout == 1
	varargout = {p};
else
	% Compute r and p_conv
	[r, p_conv] = corr(RDMa, RDMb, 'type', correlationType, 'rows', 'pairwise');
	
	if nargout == 2
		varargout = {r, p};
	elseif nargout == 3
		varargout = {r, p, p_conv};
	end%if
end%if:nargout
