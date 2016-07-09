function [sigma,shrinkage]=covdiag(x,df)
% function [sigma,shrinkage]=covdiag(x,df)
% Regularises the estimate of the covariance matrix accroding to the 
% optimal shrikage method as outline in Ledoit& Wolf (2005). 
% Shrinks towards diagonal matrix
% INPUT: 
%       x     (t*n): t observations on p random variables
%       df    (1*1): optional: degrees of freedom, if not t-1 
% OUTPUT: 
%       sigma (n*n): invertible covariance matrix estimator
%       shrinkage  : Applied shrinkage factor 
% 

% de-mean returns
[t,n]=size(x);
x=bsxfun(@minus,x,sum(x)./t);  % Subtract mean of each time series fast  

% compute sample covariance matrix
if nargin<2 || isempty(df) 
    df=t-1; 
end; 

sample=(1/df).*(x'*x);

% compute prior
prior=diag(diag(sample));

% compute shrinkage parameters
d=1/n*norm(sample-prior,'fro')^2;
y=x.^2;
r2=1/n/df^2*sum(sum(y'*y))-1/n/df*sum(sum(sample.^2));

% compute the estimator
shrinkage=max(0,min(1,r2/d));
sigma=shrinkage*prior+(1-shrinkage)*sample;

