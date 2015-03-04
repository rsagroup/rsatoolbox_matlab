%  splitStringAtCharacter is a function which takes a single string and a char in, and returns a cell array of the sub-strings between the delimited character.  If the delimiter is not present in the string, a char with just the original string is returned.
%
%  Cai Wingfield 1-2010

function splitOut = splitStringAtCharacter(stringIn, delimiter)

splitIndices = strfind(stringIn, delimiter);

if isempty(splitIndices) % If there are no splits...
	splitOut = {stringIn}; % Return the whole thing in a cell
else
	numberOfSplits = numel(splitIndices);
	splitOut = {};
	for splitPoint = 1:numel(splitIndices);
		splitOut = [splitOut; {stringIn(1:splitIndices(splitPoint)-1)}];
		lengthBeforeDeletion = length(stringIn);
		stringIn = stringIn(splitIndices(splitPoint)+1:end);
		lengthAfterDeletion = length(stringIn);
		deletedLength = lengthBeforeDeletion - lengthAfterDeletion;
		splitIndices = splitIndices - deletedLength;
	end%for
	splitOut = [splitOut; {stringIn}];
end%if

end%function
