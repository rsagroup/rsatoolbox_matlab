function M = permuteRows(M, p)
% N = permuteColumns(M, p)
% M is a matrix, p is a permutation
%
% CW 5-2010

%% Setup local constants
[nR nC] = size(M);
% nR = nRows
% nC = nColumns

%% Check for errors
if max(size(p)) ~= numel(p)
	error('permuteColumns:permutationVector', 'The permutation must be a vector.');
elseif nR ~= length(p)
	error('permuteColumns:wrongSizes', 'The length of the permutation must be the same as the number of rows of the matrix.');
end%if

%% Do the reorder
% iM = indicesMatrix
iM = repmat(p, 1, nC) + kron((0:nC-1)*nR, ones(1, nR));
M(:) = M(iM);
