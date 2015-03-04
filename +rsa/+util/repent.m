function M = repent(M,m,n)
% M = repent(M,m,n)
%
% repent is a function like repmat; but instead of tiling a matrix m x n, it
% tiles each entry m x n.
%
% Example
%
% A = [1 2; 3 4];
% repmat(A, 2, 3)
% 	1	2	1	2	1	2
%	3	4	3	4	3	4
%	1	2	1	2	1	2
%	3	4	3	4	3	4
% repent(A, 2, 3)
%	1	1	1	2	2	2
%	1	1	1	2	2	2
%	3	3	3	4	4	4
%	3	3	3	4	4	4
%
% (did I get those dimensions the right way round?)
%
% CW 5-2010

f = ones(m,n); % f = frame
M = kron(M, f);
