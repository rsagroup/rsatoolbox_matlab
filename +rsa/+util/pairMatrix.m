function C=pairMatrix(K)
% function C=pairMatrix(K)
% function to make the pairs 
% C= indicatorMatrix('allpairs',[1:K]); 
% 
% Joern Diedrichsen 2016
C=zeros(K*(K-1)/2,K);
where = 1; 
for i=1:K
    C(where:where+K-i-1,i)=1; 
    C(where:where+K-i-1,i+1:end)=-eye(K-i); 
    where = where+K-i; 
end; 
