function [Y,smoothedYfilename]=spatiallySmooth4DfMRI_mm(Y_OR_Yfilename,volSize_OR_mask,gaussianKernelFWHM_mm,voxSize_mm)
% spatially smooths a 4D fMRI data set stored in variable Y in matlab data file
% Yfilename.mat as a time-by-voxel matrix of dimensions specified in triple
% volSize = [xwidth, ywidth, zwidth].
%
% Edit by CW 6-2010: supress output to command line and graphical progress bar.
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

monitor = false;

gaussianKernelFWHM_vx=gaussianKernelFWHM_mm./voxSize_mm;
gaussianKernelSD_vx=gaussianKernelFWHM_vx/sqrt(8*log(2));

if ischar(Y_OR_Yfilename)
    load(Y_OR_Yfilename,'Y');
else
    Y=Y_OR_Yfilename;
    clear Y_OR_Yfilename;
end

if numel(volSize_OR_mask)==3
    mask=true(volSize_OR_mask);
else
    mask=volSize_OR_mask;
end

nTime=size(Y,1);


boxWidth_vx=ceil(3*gaussianKernelSD_vx)*2+1;
% 2 sds: exp(-(2)^2/2) = 0.1353
% 3 sds: exp(-(3)^2/2) = 0.0111, enough!

if gaussianKernelFWHM_mm>0
    if monitor, h=waitbar(0,['Smoothing current data set with Gaussian Kernel of ',num2str(gaussianKernelFWHM_mm),'-mm FWHM...']); end
    
    tic;
    for t=1:nTime
       if monitor, waitbar(t/nTime,h); end
        cSpatialVol=vec2map(Y(t,:),mask); %current spatial volume
        %cSpatialVol=reshape(Y(t,:),volSize); %current spatial volume

        cSpatialVol=smooth3(cSpatialVol,'gaussian',[boxWidth_vx(1) 1 1], gaussianKernelSD_vx(1));
        cSpatialVol=smooth3(cSpatialVol,'gaussian',[1 boxWidth_vx(2) 1], gaussianKernelSD_vx(2));
        cSpatialVol=smooth3(cSpatialVol,'gaussian',[1 1 boxWidth_vx(3)], gaussianKernelSD_vx(3));
        % cSpatialVol=smooth3(cSpatialVol,'gaussian',boxWidth_vx boxWidth_vx boxWidth_vx, [gaussianKernelSD_vx]); %SLOW EQUIVALENT

        Y(t,:)=map2vec(cSpatialVol,mask);
        %Y(t,:)=cSpatialVol(:)';
    end
    if monitor, close(h); end
    if monitor, disp(['Time required for smoothing data set with Gaussian Kernel of ',num2str(gaussianKernelFWHM_mm),'-mm FWHM:']); end
    t = toc;
end

if exist('Y_OR_Yfilename','var')
    smoothedYfilename=[Y_OR_Yfilename,'_gaussFWHM',num2str(gaussianKernelFWHM_vx),'vox'];
    save([smoothedYfilename,'.mat'],'Y','mask');
else
    smoothedYfilename=[];
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE AS STC FILES
sx=volSize(1); sy=volSize(2); sz=volSize(3);

% what follows is equivalent to executing save3DVoxelTimeMatAsSTCs(VT, 128, 128, simVolSizeZ, 'simVol-');
% but avoids memory problems due to passing the huge matrix VT.

% save3DVoxelTimeMatAsSTCs(VT, 128, 128, simVolSizeZ, 'simVol-');
% saves the voxel time matrix VT with 3 spatial dimensions as a set of brainvoyager STC file.

% sx and sy are the horizontal and vertical sizes of the slice, respectively.
% sz is the number of slices.
% the vertical size of VT must equal sx*sy*sz.
% sz is the slowest and sx the fastest varying dimension in the vectorized representation of space.

maximum=max(max(Y))
minimum=min(min(Y))
newRange=400;
newMinimum=100;

sliceSize=sx*sy;

for z=1:sz
    fwriteid = fopen([smoothedYfilename,'-',int2str(z),'.stc'],'w');
    count = fwrite(fwriteid,sx,'int16');
    count = fwrite(fwriteid,sy,'int16');
   
    for t=1:size(Y,1)
        YsingleSliceAsShort=(Y(t,1+(z-1)*sliceSize:z*sliceSize)-minimum)/(maximum-minimum)*newRange+newMinimum;
        count = fwrite(fwriteid,YsingleSliceAsShort,'int16');
    end
    
    status = fclose(fwriteid);
end

end%function
