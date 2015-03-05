function M = permuteColumns(M, p)
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
elseif nC ~= length(p)
	error('permuteColumns:wrongSizes', 'The length of the permutation must be the same as the number of columns of the matrix.');
end%if

%% Do the reorder
% iM = indicesMatrix
iM = repmat(1:nR, 1, nC) + kron((p-1)*nR, ones(1, nR));
M(:) = M(iM(:));

end%function
