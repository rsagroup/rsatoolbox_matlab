function T = readMask(mFiles)
% function rsa_readMask(mDir)
% Read the functional mask images for the subject.
%
% INPUTS
%   mFiles: full file path for the functional mask
%           Optionally, an Nx1 cell array of file-names can be provided to
%           load N functional mask files (requires full path)
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
% 10/2017- Added support for loading a single mask file

import rsa.spm.*

if iscellstr(mFiles)
    % if a whole cell array of file names is provided
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
else
    H       = spm_vol(mFiles);
    mask    = spm_read_vols(H);

    T.dim   = H.dim;
    T.mat   = H.mat;
    T.mask  = mask;
end;

