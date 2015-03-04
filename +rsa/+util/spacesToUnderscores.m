%  spacesToUnderscores is a function based on Niko Kriegeskorte's deunderscore
%  function.  It takes an incomming string and replaces all spaces (which can't
%  be in dynamic field names) with underscores (which can)
%
%  Cai Wingfield 2-2010

function stringORstringInCell=spacesToUnderscores(stringORstringInCell)

if iscell(stringORstringInCell)
    for lineI=1:numel(stringORstringInCell)
        line=stringORstringInCell{lineI};
        line(line==32)='_';
        stringORstringInCell{lineI}=line;
    end
else
    stringORstringInCell(stringORstringInCell==32)='_';
end

end%function
