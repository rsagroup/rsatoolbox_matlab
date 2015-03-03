function b=brightness(RGBrows)
% given triples of RGB values, returns the overall brightness vector. This
% is performed on each row of the input (RGBrows) independently.
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

RGBweights=[.241 .691 .068]';

b=sqrt(RGBrows*RGBweights);

