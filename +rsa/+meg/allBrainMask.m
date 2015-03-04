% This function generates dummy index masks fro the case of running whole
% brain analysis (a special mask with all the vertices). Aim of the
% function is to keep homogeneity across all functions and do masking in
% more standard way.
%
% [indexMasks =] allBrainMask(userOptions)
%
% IZ 3/12

function indexMasks = allBrainMask(userOptions)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

for chirality = 1:2
        switch chirality
            case 1
                chi = 'L';
            case 2
                chi = 'R';
        end%switch:chirality
        maskName = ['all_' lower(chi) 'h'];
		indexMasks.(maskName).maskIndices = 1:userOptions.targetResolution; % update IZ 02/12 previously: vertexMask + 1;
		indexMasks.(maskName).timeIndices = userOptions.temporalSearchlightLimits;
        indexMasks.(maskName).chirality = chi;
end

cd(returnHere); % And go back to where you started
