function string=compactPvalueString(p)
% given a numeric p-value returns a string that expresses it in scientific format.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

exponent=floor(log10(p));
string=[num2str(round(p*10^(-exponent+1))/10),'e',num2str(exponent)];