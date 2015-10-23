function L = defineSearchlight(Structures,Mask,varargin)
% function L = rsa.defineSearchlight(S,M)
%
% INPUTS
%   Structures:     Nx1 cell array containing the structures over which to
%                   estimate the the searchlights. Structures can either be
%                   volumetric regions of interest or surfaces. 
%
%       Volumetric regions of interest must contain the following fields:
%       - dim       1x3 vector with volume dimensions (in voxels)
%       - mat       4x4 affine transformation matrix that aligns mask with
%                   anatomical image
%       - data      vector or array with PROD(VOLDEF.voxsize)
%                   elements (logical or numerical) with brain mask
%
%       Surfaces must contain the following fields:
%       - white     vertices belonging to the white surface
%       - pial      vertices belonging to the corresponding pial surface
%       - topo      topology for the white/pial surface
%
%   Mask:       A functional mask struction for with the fields
%       - dim       1x3 vector with volume dimensions (in voxels)
%       - mat       4x4 affine transformation matrix that aligns mask with
%                   anatomical image
%       - data      vector or array with PROD(VOLDEF.voxsize)
%                   elements (logical or numerical) with brain mask
%
% OPTIONS/VARARGIN: 
%   sphere:         definition of searchlight sphere, either one of:
%                   [R]   scalar R means fixed radius R (in mm)
%                   [R C] use maximum radius R to find approximately C voxels
%   linedef:        1x3 vector which defines
%                       - linedef(1) is the number of steps between pial and white surfaces
%                       - linedef(2) is the starting distance from the white surface
%                       - linedef(3) is the ending distance from the pial surface
%   distance:       metric used to calculate distances of voxels from
%                   center voxel
%   progressstep    number of iterations after which progress is reported
%
% OUTPUT:
%       LI:         nVox x 1 cell array with linear voxel indices
%       voxmin:     nVox x 3 Minimal voxel coordinate in x,y,z direction
%       voxmax:     nVox x 3 maximal voxel coordinate in x,y,z direction
%       voxel:      nVox x 3 Matrix of I,J,K voxel coordinates or 
%                   nVox x 1 vector of linear voxel indices for the centers of search-lights 
%       structure:  Defines which of the structures the searchlight voxels
%                   belong to

%% 1. Input checking
if isempty(Structures) && isempty(Mask)
    error('provide functional mask image as a minimum input');    
end;
if isempty(Mask) && length(Structures)==1 && isfield(Structures{1},'mat')
    % run a simple volumetric searchlight
    L = defineSearchlight_volume(Structures{1},[],varargin{:});
    L.structure = ones(size(L.LI,1),1);
else
    %% 2. Mark volumetric and surface based structures
    nStruct     = length(Structures);
    whichStruct = zeros(nStruct,1);
    for i=1:nStruct
        if isfield(Structures{i},'mat')
            whichStruct(i) = 1;
        else 
            whichStruct(i) = 2;
        end;
    end;


    %% 3. Loop over the volumetric searchlights first
    idxVol = find(whichStruct==1);
    L = [];
    exclMask = Mask;
    for i=1:length(idxVol)
        [Li,exclMask]   = defineSearchlight_volume(Structures{idxVol(i)},exclMask,varargin{:});
        Li.structure    = i*ones(size(Li.LI,1),1);
        L = addstruct(L,Li);
    end;


    %% 4. Loop over the surface searchlights
    if isempty(L)
        exclMask = Mask;
        maxL = 0;
    else 
        maxL = max(L.structure);
    end;

    idxVol = find(whichStruct==2);
    for i=1:length(idxVol)
        exclMask.mask   = exclMask.data;
        [Li,exclMask]   = defineSearchlight_surface(Structures{idxVol(i)},exclMask,varargin{:});
        Li.structure    = (maxL+i)*ones(size(Li.LI,1),1);
        L = addstruct(L,Li);
    end;
end;