function varargout = replaceWildcards(stringIn, varargin)
%
% stringOut = replaceWildcards(stringIn[, wildcard1, replacement1, wildcard2, replacement2, ...])
%
% replaceWildcards is a function which takes a string "stringIn" and
% (recursively) replaces all occurrences of {wildcard1, wildcard2, ...} with
% {replacement1, replacement2, ...} until no instances of any wildcards remain.
%
% WARNING: Doing something stupid like setting replacement1 = wildcard1 or any
%          other hilarious recursive replacement set WILL cause this function to
%          recurse until either MATLAB kills it or your computer crashes.  You
%          have been warned!
%
% Cai Wingfield 5-2010
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

%% Deal with inputs

if mod(numel(varargin), 2) == 0
	nWildCards = numel(varargin)/2;
	for i = 1:nWildCards
		wildCard{i} = varargin{2*i-1};
		mildCard{i} = varargin{2*i};
	end%for:i
else
	error('replaceWildcards:wrongNumberOfArguments', ['Wrong number of arguments [' num2str(nargin) '] (must be odd).']);
end%if:odd nargin

stringOut = stringIn;

replacementsWereMade = false;

% Make replacements once
for i = 1:nWildCards
	stringOut2 = strrep(stringOut, wildCard{i}, mildCard{i});
	if ~strcmp(stringOut, stringOut2)
		replacementsWereMade = true;
	end%if
	stringOut = stringOut2;
	clear stringOut2;
end%for:i

% If replacements where made, recursively keep doing it until they're not.
while replacementsWereMade
	replacements = false(1, nWildCards);
	for i = 1:nWildCards
	[stringOut replacements(i)] = replaceWildcards(stringOut, wildCard{i}, mildCard{i});
	end%for:i
	replacementsWereMade = any(replacements);
end%while:replacementsWereMade

if nargout == 1
	varargout{1} = stringOut;
elseif nargout == 2
	varargout{1} = stringOut;
	varargout{2} = replacementsWereMade;
end

end%function
