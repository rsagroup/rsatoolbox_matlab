function L = defineSearchlight(Structures,Mask,varargin)
% function L = rsa.defineSearchlight(S,M)
%
% INPUTS
%   Structures:     Nx1 cell array containing the structures over which to
%                   estimate the searchlights. Structures can either be
%                   volumetric regions or surfaces. Voxel contained in these
%                   regions will serve as centers for the searchlights.
%
%     Volumetric regions of interest must contain the following fields:
%       - dim       1x3 vector with volume dimensions (in voxels)
%       - mat       4x4 affine transformation matrix that aligns mask with
%                   anatomical image
%       - data      vector or array with PROD(VOLDEF.voxsize)
%                   elements (logical or numerical) with brain mask
%       - sphere    definition of searchlight sphere, either one of:
%                   [R]   scalar R means fixed radius R (in mm)
%                   [R C] use maximum radius R to find approximately C voxels
%     Surfaces must contain the following fields:
%       - white     vertices belonging to the white surface
%       - pial      vertices belonging to the corresponding pial surface
%       - topo      topology for the white/pial surface
%       - sphere    definition of searchlight circle, either one of:
%                   [R]   scalar R means fixed radius R (in mm)
%                   [R C] use maximum radius R to find approximately C voxels
%       - linedef   1x3 vector which defines the cortical depth at which we sample
%                   Voxels for the search light.
%                       - linedef(1) is the number of steps between pial and white surfaces
%                       - linedef(2) is the lower surface (0 = white surface)
%                       - linedef(3) is the upper surface (1 = for pial surface)
%                   Thus linedef [0 0.5 0.5] samples a single mid-layer
%                        linedef [4 0 1.3] samples from the white to 30% outside the pial layer
%       - distance  metric used to calculate distances of voxels from
%                   center voxel: "" or "geodesic" (default)
%
%   Mask:       A functional mask define the voxels that have data that can become part of the searchlight.
%               For Volume-based search lights these are assumed to be the same.
%       - dim       1x3 vector with volume dimensions (in voxels)
%       - mat       4x4 affine transformation matrix that aligns mask with
%                   anatomical image
%       - mask      vector or array with PROD(VOLDEF.voxsize)
%                   elements (logical or numerical) with brain mask
%
% OPTIONS/VARARGIN:
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

Opt.progressstep = 100;   % Progress step reporting 
Opt.writeMask     = 0;    % Save the full mask image for included voxels? 

rsa.getUserOptions(varargin,Opt,{'progressstep','saveMask'}); 

%% 1. Input checking
if isempty(Structures) && isempty(Mask)
    error('provide functional mask image as a minimum input');
end;

% If it is a Single Structure - wrap as cell array
if (isstruct(Structures))
    Structures={Structures};
end;

% If it's a single searchlight and no mask is given - just make the
% volume-based searchlight the mask.
if isempty(Mask) && length(Structures)==1 && isfield(Structures{1},'mat')
    Mask = Structures{1};
end;

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
numStruct = 1; 
for i=1:length(idxVol)
    currentOpt=Opt; 
    if (isfield(Structures{idxVol(i)},'sphere')); 
        currentOpt.sphere=Structures{idxVol(i)}.sphere;
    end; 
    [Li,exclMask]   = rsa.defineSearchlight_volume(Structures{idxVol(i)},exclMask,currentOpt);
    Li.structure    = numStruct*ones(size(Li.LI,1),1);
    L = addstruct(L,Li);
    numStruct = numStruct+1; 
end;


%% 4. Loop over the surface searchlights
idxVol = find(whichStruct==2);
for i=1:length(idxVol)
    currentOpt=Opt; 
    if (isfield(Structures{idxVol(i)},'sphere')); 
        currentOpt.sphere=Structures{idxVol(i)}.sphere;
    end; 
    if (isfield(Structures{idxVol(i)},'linedef')); 
        currentOpt.sphere=Structures{idxVol(i)}.linedef;
    end; 
    if (isfield(Structures{idxVol(i)},'distance')); 
        currentOpt.sphere=Structures{idxVol(i)}.distance;
    end; 
    [Li,exclMask]   = rsa.defineSearchlight_surface(Structures{idxVol(i)},exclMask,currentOpt);
    Li.structure    = (numStruct)*ones(size(Li.LI,1),1);
    L = addstruct(L,Li);
    numStruct = numStruct+1; 
end;
