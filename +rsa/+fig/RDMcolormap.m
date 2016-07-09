function cols=RDMcolormap
% this function provides a convenient colormap for visualizing
% dissimilarity matrices. it goes from blue to yellow and has grey for
% intermediate values.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

nCols = 256;
%% blue-cyan-gray-red-yellow with increasing V (BCGRYincV)
anchorCols=[0 0 1
            0 1 1
            .5 .5 .5 
            1 0 0
            1 1 0];

anchorCols_hsv=rgb2hsv(anchorCols);
incVweight=1;
anchorCols_hsv(:,3)=(1-incVweight)*anchorCols_hsv(:,3)+incVweight*linspace(0.5,1,size(anchorCols,1))';

brightness(anchorCols);
anchorCols=hsv2rgb(anchorCols_hsv);

cols=colorScale(anchorCols,nCols);
cols1=cols;

% figure(1); colormap(cols);
% cm=colormap; hsvcm=rgb2hsv(cm);
% subplot(3,1,3); cla; plot(brightness(colormap),'r'); hold on; plot(hsvcm(:,3),'k'); axis tight; legend({'brightness','V'});
% th=addHeading('blue-cyan-gray-red-yellow with increasing V');

end%function
