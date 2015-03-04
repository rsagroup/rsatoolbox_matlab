function [adjacents, passed] = getadjacent(str1, int, hashtab)

% by Li Su and Andy Thwaites

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

adjacentsbelow = [];
adjacents = [];
passed = [];
if int==1
   adjacents = hashtab.get(str1);
   passed = [1];
else
   [adjacentsbelow, passed] = getadjacent(str1, int-1, hashtab);
   for j = 1:length(adjacentsbelow)
      adjacents = [adjacents; hashtab.get(num2str(adjacentsbelow(j)))];
   end
   adjacents = unique(adjacents);
   passed = [passed; adjacentsbelow];
   for j = length(adjacents):-1:1
      if(any(find(passed == adjacents(j))))
          adjacents(j)=[];
      end
   end
end

adjacents = unique(adjacents);
passed = unique(passed);
