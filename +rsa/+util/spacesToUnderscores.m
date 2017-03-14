function stringORstringInCell=spacesToUnderscores(stringORstringInCell)
%
%  spacesToUnderscores is a function based on deunderscore.m
%  It takes an incoming string and replaces all spaces (which
%  as underscores.
%
%  based on underscoresToSpaces written by Cai Wingfield 11-2009
% 
%  Ian Charest 3-2017
%_________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

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