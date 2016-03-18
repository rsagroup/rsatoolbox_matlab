function [omega,R2,loglike]=fitModelOLS(Model,dist,varargin);
% function [omega,theta]=rsa_fitModelOLS(Model,Data,varargin);
% Does a linear fit to some data
% INPUT:
%   Model:     numDist x numReg:  Design matrix or model structure 
%   dist:      numObs x numDist:  Data, with obserations stacked in rows
% OUTPUT:
%   omega:     estimates of the regression coefficients, 1 row per observation 
%   R2:        R2 of the fit: If intercept is set to 0, this is as proportion of the
%              the total-sums-of-squares sum(d.^2). 
%              If intercept is set to 1, this is given in respect to sum((d-mean(d)).^2) 
%   loglike    loglikelihood of the data uner the model, using an i.i.d
%              normal approximation and a point ML-etimate for the variance
%              
import rsa.util.*;
import rsa.*;

Opt.intercept = 0;          % Implicitly remove intercept of the distances? Default is no
Opt=rsa.getUserOptions(varargin,Opt);

[numObs,numDist]=size(dist); 
numCond = ceil(sqrt(numDist*2));
if (numCond*(numCond-1)/2~=numDist) 
    error('bad size of distances'); 
end; 

% Deal with different ways of specifying the model 
if (isstruct(Model))
    if (length(Model)>1)        % Structure array 
        X=vertcat(Model(:).RDM)'; 
    else 
        X=permute(Model.RDM,[2 1 3]);
    end;
else 
    X=Model; 
end; 
        
numReg = size(X,2); 

% Get the data 
Y=dist'; 

% If intercept: remove from data and design matrix 
% the R2 fit measure will then also reflect this 
if Opt.intercept 
    X = bsxfun(@minus,X,mean(X,1)); 
    Y = bsxfun(@minus,Y,mean(Y,1)); 
end; 

% Do the regression analysis 
omega =(X'*X)\(X'*Y);       
RSS   = sum((Y-X*omega).^2,1);
R2    = 1-RSS./sum(Y.^2,1);

if (nargout>2) 
    sigma = RSS / numDist;                                                            % Maximum liklihood
    loglike  = -numDist/2*log(2*pi)-numDist/2*log(sigma)-1/(2*numDist); % Log likelihood of the data at maximum likelihood 
    loglike  = loglike'; 
end; 
% Transpose return arguments 
omega = omega'; 
R2    = R2'; 