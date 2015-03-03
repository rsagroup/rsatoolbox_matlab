function v = randomlyPermute(v)
% randomly permutes the elements of the input vector (v).
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council
v = v(randomPermutation(numel(v)));