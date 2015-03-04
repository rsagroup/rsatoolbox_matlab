%  underscoresToSpaces is a function based on Niko Kriegeskorte's deunderscore
%  function.  It takes an incomming string and replaces all dashes (which
%  can't be used as field names) as underscores (which can).
%
%  Li Su 2-2012

function stringORstringInCell=dashToUnderscores(stringORstringInCell)

if iscell(stringORstringInCell)
    for lineI=1:numel(stringORstringInCell)
        line=stringORstringInCell{lineI};
        line(line==45)='_';
        stringORstringInCell{lineI}=line;
    end
else
    stringORstringInCell(stringORstringInCell==45)='_';
end

end%function
