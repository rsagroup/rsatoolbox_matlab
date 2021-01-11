function [L,exclMask] = defineSearchlight_volume(ROIMask,Mask,varargin)
% function L = rsa.defineSearchlight(Mask,Opt)
% Defines a volumetric searchlight using all voxels that are > 0 in both 
% ROIMask and the Mask as centres for the searchlights and as potential voxels 
% to be included in the searchlight. 
% The Mask also defines the voxel space in which the searchlight will be
% defined. 
%
% INPUTS
%   ROIMask:    A mask struction for the region of interest with the fields
%       - dim       1x3 vector with volume dimensions (in voxels)
%       - mat       4x4 affine transformation matrix that aligns mask with
%                   anatomical image
%       - data      vector or array with PROD(VOLDEF.voxsize)
%                   elements (logical or numerical) with brain mask
%
%   Mask:       A functional mask structure that defineds the voxel
%               available  for the searchlight and the voxel space. This
%               will be usually the mask from the first-level GLM. 
%       - dim       1x3 vector with volume dimensions (in voxels)
%       - mat       4x4 affine transformation matrix that aligns mask with
%                   anatomical image
%       - mask      vector or array with PROD(VOLDEF.voxsize)
%                   elements (logical or numerical) with brain mask
%
% OPTIONS/VARARGIN: 
%   	sphere:     definition of searchlight sphere, either one of:
%                   [R]   scalar R means fixed radius R (in mm)
%                   [R C] use maximum radius R to find approximately C voxels
%
% OUTPUT:
%       Searchlight structure with the fields: 
%          LI:      nVox x 1 cell array with linear voxel indices
%          voxmin:  nVox x 3 Minimal voxel coordinate in x,y,z direction
%          voxmax:  nVox x 3 maximal voxel coordinate in x,y,z direction
%          voxel:   nVox x 3 Matrix of I,J,K voxel coordinates or 
%                   nVox x 1 vector of linear voxel indices for the centers of search-lights 
%
% EXAMPLE:
%   Define a volumetric searchlight over the volume of the cerebellum, as
%   defined by a mask in anatomical space, using all the voxels available
%   in the functional mask. Each searchlight will have 100 voxels 
%   ROI         = spm_vol('anatomical/p03_anatomical_cerebellum.nii');
%   ROI.data    = spm_read_vols(ROI);
%   M           = rsa.readMask('glm/p03/mask.img');
%   L           = rsa.defineSearchlight_volume(ROI,M,'sphere',[50 100]);
%
% 2/2015 - Joern Diedrichsen & Naveed Ejaz 

import rsa.spm.*

%% 1. Checking inputs
Opt.sphere    = [30 160]; % 30mm radius volumetric sphere, 160 voxels
Opt.writeMask = 0;        % Optional writing of the mask.  
Opt         = rsa.getUserOptions(varargin,Opt); 

if isempty(ROIMask) && isempty(Mask)
    error('provide functional mask image');
end;

if isempty(Mask)
    Mask = ROIMask;
elseif isempty(ROIMask)
    ROIMask = Mask;
end;

if ~isequal(size(ROIMask.mat),[4 4])
        error('affine matrix for mask %i should be 4x4',n);
end;
if ~isequal(size(ROIMask.data),ROIMask.dim)
        error('affine matrix for mask %i should be 4x4',n);
end;


%% 2. Getting sphere definition, and setting up reference mask
radius      = Opt.sphere(1);
fixedradius = numel(Opt.sphere)==1;
if ~fixedradius
    targetvoxelcount = Opt.sphere(2);
end;


%% 3. Estimate searchlight voxel indices over the provided region
% - estimate valid voxels given region and functional masks
Vin(1) = Mask; 
Vin(2) = ROIMask; 
Vo     = Vin(1); 
Vo.fname =  'inclMask.nii'; 
% inclMask        = spm_imcalc({Mask.fname; ROIMask.fname},'inclMask.nii','i1.*i2');
inclMask        = spm_imcalc(Vin,Vo,'i1.*i2');
inclMask.data   = spm_read_vols(inclMask);

%	- calculate voxels of interest and center voxels
mask        = logical(inclMask.data(:)>0);
centeridxs  = find(mask);   
c_centers   = surfing_inds2subs(inclMask.dim,centeridxs)'; 
c_voxels    = surfing_inds2subs(inclMask.dim,find(mask))'; 
centers     = uint32(c_centers);
voxels      = uint32(c_voxels);
ncent       = size(centers,2); 

%   - preallocate the matrices/cell arrays
li          = cell(ncent,1);        % linear indices for voxels
n           = zeros(ncent,1);       % Number of voxels 
rs          = zeros(ncent,1);       % Searchlight radius 
voxmin      = zeros(ncent,1);       % bottom left voxel
voxmax      = zeros(ncent,1);       % top right voxel


%% 4. Estimate linear indices, voxmin/voxmax for a sphere centered at each center index
spm_progress_bar('Init',100);
for k=1:ncent
    ds = surfing_eucldist(c_centers(:,k),c_voxels);
    if fixedradius
        a       = voxels(:,ds<radius);
        rs(k,1) = radius; 
    else 
        i       = find(ds<radius(1));
        [dss,j] = sort(ds(i));
        indx    = min(targetvoxelcount,length(i)); % In case there are not enough voxels within maximal radius
        a       = voxels(:,i(j(1:indx))); 
        rs(k,1) = dss(indx); 
    end; 
    n(k,1)          = size(a,2);
    voxmin(k,1:3)   = min(a,[],2)';
    voxmax(k,1:3)   = max(a,[],2)';
    li{k,1}         = surfing_subs2inds(inclMask.dim,a(1:3,:)')';
    n(k,1)          = numel(li{k});

    spm_progress_bar('set',(k/ncent)*100);
end;


%% 5. Setting output
L.LI        = li;
L.voxmin    = voxmin;
L.voxmax    = voxmax;
L.voxel     = centeridxs;


%% 6. Cleanung & writing out exclusion mask
Vin(1) = Mask; 
Vin(2) = inclMask; 
if (Opt.writeMask) 
    Vo     = Vin(1); 
    Vo.fname =  'exclMask.nii'; 
    exclMask        = spm_imcalc(Vin,Vo,'((i1>0)-(i2>0)>0)');
    exclMask.data   = spm_read_vols(exclMask);
end; 
try
    delete('inclMask.nii');
catch
end; 