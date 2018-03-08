function stringORstringInCell=deunderscore(stringORstringInCell)
% gets an input string and replaces the underscores ('_') with hyphens
% ('-').This avoids the problems with displaying the names in which the
% first character after the underscore would be displayed as an index. 
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

% replace underscores
if iscell(stringORstringInCell)
    for lineI=1:numel(stringORstringInCell)
        line=stringORstringInCell{lineI};
        line(line==95)='-';
        stringORstringInCell{lineI}=line;
    end
else
    stringORstringInCell(stringORstringInCell==95)='-';
end

end%function
