function y=surfing_uniqueidxsperrow(x)
% returns unique indices per row
%
% UNIQUEIDXS=SURFING_UNIQUEIDXSPERROW(IDXS)
% 
% INPUT:
%   IDXS:       NxP matrix, for N rows with P indices each.
% OUTPUT:
%   UNIQUEIDXS: NxQ matrix, Q<=P, where each row contains the unique
%               values in the corresponding row from IDXS. 
%               Only positive values in IDXS are returned, and
%               each row may contain zeros in case duplicate values
%               were found in that row.
% 
% This function has also been implemented in C, which runs a lot faster 
% (output is slightly different but conforms to the specification).
% If MEX is set up, run "mex surfing_uniqueidxsperrow.c" for increased
% speed.
%
% NNO May 2010

[nrows,ncols]=size(x);

y=zeros(nrows,ncols);
for k=1:nrows
    unqxk=unique(x(k,x(k,:)>0)); %only consider positive elements, and take the unique values
    y(k,1:numel(unqxk))=unqxk;   %store them in y
end

% the number of necessary columns
yncols=find(sum(y>0)==0,1)-1;

% if all columns have a non-zero element then yncols is empty, and y is kept as is.
if ~isempty(yncols) 
    y=y(:,1:yncols);
end

