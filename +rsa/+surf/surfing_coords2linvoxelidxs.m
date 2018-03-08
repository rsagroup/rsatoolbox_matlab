function linidxsrs=surfing_coords2linvoxelidxs(coords,voldef,dim,mask)
% maps coordinates to linear voxel indices
%
% I=SURFING_COORDS2VOXELIDXS(C,VOLDEF)
% INPUTS: 
%   C:       A 3xN matrix or 3xPxQ array with (x,y,z) coordinates 
%   VOLDEF:  Struct with fields .mat (4x4 voxel to coordinate transformation matrix 
%            rom the images to be sampled (1-based)) and .dim (1x3
%            volume dimensions in voxels)
%            Optionally a field .mask may be specified (Nx1 or IxJxK) that 
%            specifies which voxels should be included. 
% OUTPUT:
%   I:       Nx1 vector or PxQ matrix with linear voxel indices
%            Indices outsize the volume (and if voldef.mask is specified, 
%            unselected voxels) are set to NaN.
% 
% * Alternative input is SURFING_COORDS2VOXELIDXS(C,MAT,DIM[,MASK])
% * To obtain unique indices, run SURFING_UNIQUEIDXSPERROW
% 
% TW,NNO May 2010
%
% See also SURFING_NODEIDXS2COORDS, SURFING_UNIQUEIDXSPERROW,
% SURFING_INDS2SUBS

% check input
if nargin<3
    if ~isstruct(voldef) || ~isfield(voldef,'mat') || ~isfield(voldef,'dim')
        error('voldef should be a struct with fields .mat and .dim');
    end
    mat=voldef.mat;
    dim=voldef.dim;
	
    
    % NNO Sep 2011 - allow for separate centermask that is different from
    % volume mask
    mask    = NaN;
    usemask = 0;
    
    if isfield(voldef,'centermask');
        mask    = voldef.centermask;
        usemask = 1;
        fprintf('using centermask\n');
    elseif isfield(voldef,'mask'); 	% see if voxel mask is specified
        mask    = voldef.mask;
        usemask = 1;
        fprintf('using mask\n');
    end

else
    mat=voldef;
	usemask=nargin>=4; 	% see if voxel mask is specified
end

if usemask
	% check mask for correct number of elements
	if numel(mask) ~= prod(dim)
		error('Mask size does not match voxel dimension');
    end
		
    if isnumeric(mask)
        mask=abs(mask)>0; % exclude values of 0 and NaN
    end
end

if ~isequal(size(mat),[4 4]), error('matrix should be 4x4'); end
if ~isequal(size(dim),[1 3]), error('dim should be 1x3 vector'); end

rs=size(coords);
if rs(1) ~= 3, error('First dimensions of coords should be 3'); end

% see if coords is 3xN, or 3xPxQ
switch numel(rs)
    case 2
        % coords is 3xN
        ncoordspernode=1;
        nverts=rs(2);
    case 3
        % coords is 3xPxQ
        ncoordspernode=rs(2); % P
        nverts=rs(3);         % Q
    otherwise
        error('coords has %d dimensions, not supported', numel(rs));
end

% map to 3xP matrix (P coordinates)
coords= reshape(coords, 3,[]);
coords=[coords;ones(1,size(coords,2))];                  

ijk=(mat\coords);  % apply transformation
ijk=round(ijk(1:3,:)); %keep (x,y,z) dimensions, and round the result to nearest integer

alllinidxs=surfing_subs2inds(dim,ijk'); % make it linear indices 
linidxs=alllinidxs; 

if usemask
    insidevolmask=~isnan(alllinidxs); % voxels inside the volume
    voxelmask=false(size(alllinidxs)); 
    voxelmask(insidevolmask)=mask(alllinidxs(insidevolmask)); % conjunction of inside the volume and selected with voxel mask
    linidxs(~voxelmask)=NaN; % set values outside the conjunction mask to NaN
end	


% % NEW
% epivoxel= find(~isnan(alllinidxs));
% goodvoxel= voldef.mask(alllinidxs(epivoxel));
% linidxs= NaN(size(alllinidxs));
% linidxs(epivoxel(find(goodvoxel>=1)))= alllinidxs(epivoxel(find(goodvoxel>=1)));
% % %END NEW


linidxsrs=reshape(linidxs,ncoordspernode,nverts)'; % reshape, make it NVERTS x NCOORDSPERNODE

