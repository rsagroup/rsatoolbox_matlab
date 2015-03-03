function addImageSequenceToAxes(ax,il)
% USAGE
%         addImageSequenceToAxes(ax,imagelabels)
% FUNCTION
%         places a set of images at equal distances along a line in a 2d
%         coordinate system, which is then mapped to the axes of the RDM
%         in ax. imagelabels is a struct with the following fields:
%            images : struct containing
%                image : RGB image matrix (x,y,3)
%                alpha: optional alpha layer for image
%            sequence : vector of custom indices for images. Optional.
%            nRows : number of rows/columns to plot labels in. Default 2.
%            transparentCol : Remove background based on rgb colour. Default undefined.
%            blackdisks : Logical. Place black disks under images. Default false.

% Default parameters
% In principle the only default we can't set is il.images
% transparentCol is left out here
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

defaultil = struct('sequence',1:length(il.images),'nRows',2, ...
	'blackdisks',false);

for fn = fieldnames(defaultil)'
	if ~isfield(il,fn{1})
		il.(fn{1}) = defaultil.(fn{1});
	end
end

nImagesInSeq=numel(il.sequence);

% Get axes
axes(ax);
da_orig = daspect;

set(ax,'layer','bottom');

xlim=get(gca,'XLim');
ylim=get(gca,'YLim');

% figure out how large the plotted RDM is
rdmlen = length(get(min(get(ax,'children')),'cdata'));

% The RDMs get plotted with a .5 offset
offset = .5;

% Offset in Z to make lines reliably appear behind faces
% Doesn't actually do anything (?)
offset_z = 500;

% Padding the axes slightly at the end helps with images clipping outside the frame
padding = 1.05;

% Add space to X and Y to make room for labels
set(ax,'xlim',[xlim(1) - il.nRows xlim(2)]);
set(ax,'ylim',[ylim(1) - il.nRows ylim(2)]);

% Get updated values
xlim=get(gca,'XLim');
ylim=get(gca,'YLim');


% Place in the space we've created
start_xy_c{1} = [xlim(1)-1 ylim(2)-offset];
end_xy_c{1} = [xlim(1)-1 ylim(1)+il.nRows+offset];

start_xy_c{2} = [xlim(2)-offset ylim(1)+1];
end_xy_c{2} = [xlim(1)+il.nRows+offset ylim(1)+1];


% Making these changes tends to mess up the axes slightly..
daspect(da_orig);

da=daspect; % need this to render image pixels square
set(gca,'XTick',[]); % switch off horizontal-axis ticks and image index numbers
hold on;

% DIMENSION
for dim = 1:2
	start_xy = start_xy_c{dim};
	end_xy = end_xy_c{dim};

	%daspect(daspect);
	% looks moot, but causes the current data aspect ratio to be preserved when the figure is resized

	%% compute image size
	% images assumed to be square
	seqVec_sqA=(end_xy-start_xy)./da(1:2); % vector in il.sequence direction (in square axis)
	imageWidth=max(abs(seqVec_sqA))/(nImagesInSeq-1)*il.nRows;
	imW_ax=imageWidth*da(1);
	imW_ay=imageWidth*da(2);

	%% prepare multi-row arrangement
	% compute vector orthogonal to il.sequence direction
	orthVec=-(null(seqVec_sqA)'*imageWidth).*da(1:2);

	% shift start and end positions to center the multiple rows on the
	% requested il.sequence line
	% start_xy=start_xy-orthVec*(il.nRows-1)/2;
	% end_xy=end_xy-orthVec*(il.nRows-1)/2;

	%% arrange the images
	for sequenceI=1:nImagesInSeq
		imageI=il.sequence(sequenceI);
		xy=start_xy+(end_xy-start_xy)/(nImagesInSeq-1)*(sequenceI-1);
		
		xy=xy+mod(sequenceI,il.nRows)*orthVec; % lateral displacement for multiple rows
		
		%[xs,ys,rgb3]=size(il.images(imageI).image); % assuming square images for now, so this isn't needed

		% If there is a transparency column, use it
		% If there is also an alpha channel, use that too
		% If there is nothing, just use whatever we have..
		st = size(il.images(imageI).image(:,:,1));

		% If there is a colour that should be made transparent
		if exist('il.transparentCol','var')
			transparent=il.images(imageI).image(:,:,1)==il.transparentCol(1) & il.images(imageI).image(:,:,2)==il.transparentCol(2) & il.images(imageI).image(:,:,3)==il.transparentCol(3);
		else
			% Once again, assuming square stimuli...
			transparent = logical(zeros(st(1),st(2)));
		end

		% If we have an alpha layer
		if isfield(il.images(imageI),'alpha')
			alpha = il.images(imageI).alpha == 0;
		else
			alpha = logical(zeros(st(1),st(2)));
		end

		% Merge colour and alpha layers
		transparentalpha = (transparent + alpha) > 0;

		% black disks underneath
		% TODO: Make these actually work
		if il.blackdisks
			angles=0:0.1:2*pi;
			X=sin(angles)*imW_ax/2+xy(1);
			Y=cos(angles)*imW_ay/2+xy(2);
			Z=-2*ones(size(X));
			patch(X,Y,Z,[0 0 0]);
		end
		
		%rectangle('Position',[xy(1)-imW_ax/2  xy(2)-imW_ay/2 imW_ax imW_ay],'Curvature',[1 1],'FaceColor','k','EdgeColor','none','EraseMode','normal');
		%plot(xy(1),xy(2),'o','MarkerSize',10,'MarkerFaceColor','k');

		% Tick line
		switch dim
			case 1
				tick_y = repmat(mean(xy(2)+[imW_ay/2, -imW_ay/2]),2,1);
				tick_x = [mean(xy(1)+[imW_ax/2, -imW_ax/2]) .5];
			case 2
				tick_x = repmat(mean(xy(1)+[imW_ax/2, -imW_ax/2]),2,1);
				tick_y = [mean(xy(2)+[imW_ay/2, -imW_ay/2]) .5];
		end

		pl(dim).x(:,sequenceI) = tick_x;
		pl(dim).y(:,sequenceI) = tick_y;
		pl(dim).z(:,sequenceI) = [offset_z offset_z];
		
		% low-level version of image function
		imh(sequenceI) = image('CData',flipdim(il.images(imageI).image,1),'XData',xy(1)+[-imW_ax/2, imW_ax/2],'YData',xy(2)+[imW_ay/2, -imW_ay/2],'AlphaData',~flipdim(transparentalpha,1),'EraseMode','normal');
		
		% high-level version of image function
		% image(xy(1)+[-imW_ax/2, imW_ax/2],xy(2)+[imW_ay/2, -imW_ay/2],il.images(imageI).image,'AlphaData',~transparentalpha,'EraseMode','normal');
		
		xlim(1)=min(xlim(1),xy(1)-imW_ax/2);
		xlim(2)=max(xlim(2),xy(1)+imW_ax/2);
		
		ylim(1)=min(ylim(1),xy(2)-imW_ay/2);
		ylim(2)=max(ylim(2),xy(2)+imW_ay/2);
		
	end

	set(gca, 'XLim',xlim, 'YLim',ylim);

end

%switch dim
	%case 1
		%pl(dim).x(:,end+1) = [offset offset];
		%pl(dim).y(:,end+1) = [start_xy(2)+offset end_xy(2)-offset];
		%pl(dim).z(:,end+1) = [offset_z offset_z];
	%case 2
		%pl(dim).x(:,end+1) = [start_xy(1)+offset end_xy(1)-offset];
		%pl(dim).y(:,end+1) = [offset offset];
		%pl(dim).z(:,end+1) = [offset_z offset_z];
%end

%start_xy = start_xy_c{dim};
% Add a main supporting line
main.x = [start_xy_c{2}(1)+offset end_xy_c{2}(1)-offset offset offset];
main.y = [offset offset start_xy_c{1}(2)+offset end_xy_c{1}(2)-offset];
main.z = [offset_z offset_z offset_z offset_z];

% Plot all the tick lines in one go
l = line([pl(1).x pl(2).x],[pl(1).y pl(2).y],[pl(1).z pl(2).z]);
set(l,'color','k','linewidth',1);
uistack(l,'bottom');

% And the main supporting line
lmain = line(main.x,main.y,main.z);
set(lmain,'color','k','linewidth',1);
uistack(lmain,'bottom');

% Add a little bit of padding
oldx = get(ax,'xlim');
oldy = get(ax,'ylim');
set(ax,'xlim',oldx * padding);
set(ax,'ylim',oldy * padding);

end%function
