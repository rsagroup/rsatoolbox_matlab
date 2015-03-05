% spawnSPMStruct
%
% Make a blank SPM metadata struct ripe for impregnation with real stuff.
%
% CW 6-2010

function V = spawnSPMStruct()

    V.fname = '';
    V.mat = eye(4);
    V.dim = [0 0 0 0];
    V.pinfo = [1;0;0];

end%function

