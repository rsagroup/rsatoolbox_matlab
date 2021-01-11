function ds=surfing_eucldist(xs,ys)
% computes the euclidian distance between points
%
% DS=SURFING_EUCLDISTANCEFROM(XS,YS) 
% INPUTS:
%   XS:   3xP matrix for P coordinates
%   YS:   3xQ matrix for Q coordinates
% OUTPUT
%   DS:   PxQ matrix, where DS(i,j) is the Euclidian distance between
%         XS(:,i) and YS(:,i)
%
% This function has also been implemented in C, which runs a lot faster.
% If MEX is set up, run "mex surfing_eucldist.c" for increased speed.
%
% NNO May 2010

[three1,n1]=size(xs);
if three1 ~= 3, error('xs should be 3xN'); end

[three2,n2]=size(ys);
if three2 ~= 3, error('ys should be 3xN'); end

% assume N1 is much less than N2; if not swap the order.
if three1>three2
    ds=surfing(eucldist(ys,xs))';
else
    ds=zeros(n1,n2);
    for k=1:n1
        % use array indexing for faster execution
        ds(k,:)=sqrt(sum((repmat(xs(:,k),1,n2)-ys).^2))';
    end
end