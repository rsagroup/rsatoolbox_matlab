%  Mex version to compute C*A*C'
%  Where C is the pairmatrix of K*(K-1)/2 x K size
%  This consitutes the main time sink in RSA-based computations,
%  for example for varianceLDC
% INPUT:
%   A: KxK matrix: single or double 
%
% OUTPUT:
%   B : K*(K-1)/2 x K matrix
%
% (c) Joern Diedrichsen