% permutations = exhaustivePermutations(n[, cap])
%
% n 			number of items to be permuted
% cap			(optional) stop after finding this many
% permutations	each row is a permutation; they'll be in lexicographic order
%
% Cai Wingfield 3-2010

function permutations = exhaustivePermutations(varargin)

	% The following algorithm generates the next permutation lexicographically after
	% a given permutation. It changes the given permutation in-place.
	% 
	%    1. Find the largest index j such that a[j] < a[j + 1]. If no such index
	%       exists, the permutation is the last permutation.
	%    2. Find the largest index l such that a[j] < a[l]. Such a l exists and
	%       satisfies j < l, since j + 1 is such an index.
	%    3. Swap a[j] with a[l].
	%    4. Reverse the sequence from a[j + 1] up to an including the final element
	%       a[n].
	%
	% After step 1, one knows that all of the elements strictly after position j
	% form a weakly decreasing sequence, so no permutation of these elements will
	% make it advance in lexicographic order; to advance one must increase a[j].
	% Step 2 finds the smallest value a[l] to replace a[j] by, and swapping them in
	% step 3 leaves the sequence after position j in weakly decreasing order.
	% Reversing this sequence in step 4 then produces its lexicographically minimal
	% permutation, and the lexicographic successor of the initial state for the
	% whole sequence.

	switch nargin
		case 1
			n = varargin{1};
			cap = inf;
		case 2
			n = varargin{1};
			cap = varargin{2};
		otherwise
			error('exhaustivePermutations:nargin', 'Only one or two arguments, please.');
	end%switch:nargin

	absoluteMaximum = 21; % Larger than this and n! is not accurately represented by a double. Does this matter?

	initialLexOrder = 1:n;

	%if n > absoluteMaximum, permutations = []; return; end%if
	if n > absoluteMaximum, warning('exhaustivePermutations:tooMany', ['21! is the largest number accurately represented by a double.\nAttempting ' num2str(n) '!, therefore, may give\nunpredictable results.']); end%if

	% The first in the list is just the initial lexical order
	permutations = initialLexOrder;
	i = 1;
	allPermutationsFound = false;

	while ~allPermutationsFound && i < cap

		currentPermutation = permutations(i,:);
		[nextPermutation, allPermutationsFound] = getNextPermutation(currentPermutation);
		permutations = [permutations; nextPermutation];
		i = i + 1;

	end%while:~allPermutationsFound
	
end%function

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%

function [nextPermutation, allPermutationsFound] = getNextPermutation(currentPermutation)

	n = numel(currentPermutation);

	allPermutationsFound = false;

	% Search through the current permutation, beginning to end, to find the largest index j such that currentPermutation(j) < currentPermutation(j+1)
	j = 0;
	for jSearch = 1:n-1
		if currentPermutation(jSearch) < currentPermutation(jSearch+1)
			j = jSearch;
		end%if
	end%for:jSearch

	% If no such j exists, the current permutation is the last permutation in lexicographical order.
	% Otherwise, find the largest index k such that currentPermutation(k) > currentPermutation(j).
	if j == 0

		nextPermutation = [];
		allPermutationsFound = true;
	
	else

		k = 0;
		for kSearch = j:n
			if currentPermutation(kSearch) > currentPermutation(j)
				k = kSearch;
			end%if
		end%for:kSearch
	
		nextPermutation = currentPermutation;
	
		% Swap currentPermutation(j) and currentPermutation(k)
	
		nextPermutation_j = nextPermutation(j);
		nextPermutation_k = nextPermutation(k);
	
		nextPermutation(j) = nextPermutation_k;
		nextPermutation(k) = nextPermutation_j;
	
		% Reverse order of currentPermutation(j+1:end)
	
		nextPermutation(j+1:end) = nextPermutation(end:-1:j+1);
	
	end%if:j == 0

end%function
