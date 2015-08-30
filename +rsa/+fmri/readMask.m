function T = readMask(mFiles)
% function rsa_readMask(mDir)
% Read the functional mask images for the subject.
%
% INPUTS
%   mFiles: Nx1 cell array with full path and file names for the subjects 
%           functional mask file
%
% OUTPUTS
%   T: Nx1 structure with the fields
%       - dim       1x3 vector with volume dimensions (in voxels)
%       - mat       4x4 affine transformation matrix that aligns mask with
%                   anatomical image
%       - mask      vector or array with PROD(VOLDEF.voxsize)
%                   elements (logical or numerical) with brain mask
%
% EXAMPLE:
%   % Read mask for subject p03
%   mDir = ~/Documents/data/rsa_surfing/glm/p03
%   T = rsa_readMask(mDir);
%
% Naveed Ejaz
% n.ejaz@ucl.ac.uk
% 2/2015

import rsa.spm.*

nSubj = length(mFiles);
for i=1:nSubj
    mFile = mFiles{i};
    if ~exist(mFile,'file')
        fprintf('Mask image %i not found ...\n',i);
    else
        % 1. Read mask image header & data
        H    = spm_vol(mFile);
        mask = spm_read_vols(H);

        % 2. Output structure
        T(i).dim  = H.dim;
        T(i).mat  = H.mat;
        T(i).mask = mask;
    end;
end;


