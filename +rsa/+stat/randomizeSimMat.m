function randomizedSimMat=randomizeSimMat(simMat)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

[v,h]=size(simMat);

if v~=h
    % convert from upper triangular form
    simMat=squareform(simMat);
end

scrambledIs=randperm(size(simMat,1));

randomizedSimMat=simMat(scrambledIs,scrambledIs);

if v~=h
    % convert to upper triangular form
    randomizedSimMat=squareform(randomizedSimMat);
end


