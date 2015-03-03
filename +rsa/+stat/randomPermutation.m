function indicesOut = randomPermutation(n)
% generates a random vector of n integer numbers from 1 to n. The numbers
% are randomly placed in the vector.
% Cai Wingfield 2-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

indicesToChoseFrom = (1:n);

for i = 1:n

	position = ceil(rand*numel(indicesToChoseFrom));
	thisIndexChoice = indicesToChoseFrom(position);
	indicesOut(i) = thisIndexChoice;
	indicesToChoseFrom = removeElement(indicesToChoseFrom, position);

end%for(i)

end%function

%% === Subfunctions =============================

function vectorOut = removeElement(vectorIn, i)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if i == 1;
	vectorOut = vectorIn(2:end);
elseif i > 1 && i < numel(vectorIn)
	vectorOut = vectorIn([1:i-1, i+1:end]);
elseif i == numel(vectorIn)
	vectorOut = vectorIn(1:end-1);
else
	error('Can''t remove an element which isn''t in the input!');
end%if

end%function
