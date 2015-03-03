function sigma=covdiag(x)

% function sigma=covdiag(x)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards diagonal matrix

% de-mean returns
[t,n]=size(x);
meanx=mean(x);
x=x-meanx(ones(t,1),:);

% compute sample covariance matrix
sample=(1/t).*(x'*x);

% compute prior
prior=diag(diag(sample));

% compute shrinkage parameters
d=1/n*norm(sample-prior,'fro')^2;
y=x.^2;
r2=1/n/t^2*sum(sum(y'*y))-1/n/t*sum(sum(sample.^2));

% compute the estimator
shrinkage=max(0,min(1,r2/d));
sigma=shrinkage*prior+(1-shrinkage)*sample;
