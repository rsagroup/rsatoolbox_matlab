function [y1,y2,y3] = rsa_affine_transform(x1,x2,x3,M)
% function [y1,y2,y3] = affine_transform(x1,x2,x3,M)
% -----------------------------------------------------------------
% Affine Transform for input stuff in any format (N-dim strcutures)
% -----------------------------------------------------------------
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);

