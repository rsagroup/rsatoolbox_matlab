function stringOut = removeSpaces(stringIn);
%
% removeSpaces will remove all spaces from a string
%
% Cai Wingfield 1-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council



if ~isempty(stringIn)
	stringOut = strrep(stringIn, ' ', '');
end
