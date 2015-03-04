function r_matrix = g_matrix(data, nConditions, nRepsPerCondition)
    % computes g_matrix for given input parameters
    % data in format (:,:,:,session)
    % nConditions: total conditions in experiment
    % nRepsPerConditions: number of time data for same condition is
    %                       repeated
    % r_matrix: correlation matrix (need to compute 1-r for Spearman corr)
    
    % based on Diedrichson et al. / NeuroImage (2011) 1665 - 1678
    % Isma Zulfiqar (2012)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

% calculating z matrix
modelBase = eye(nConditions);
z_matrix = zeros(nConditions * nRepsPerCondition, nConditions);
itr=1;
for i=1:nConditions
    for j = 1:nRepsPerCondition
        z_matrix(itr,:) = modelBase(i,:);
        itr = itr + 1;
    end
end
                    
% calculating A matrix
D.K = nConditions;
Ac={};
for i=1:D.K
    for j=i:D.K
        Ac{end+1}=zeros(D.K,D.K);
        Ac{end}(i,j)=1;
        Ac{end}(j,i)=1;
    end;
end
                    
% calculating G matrix
[G,h,u,l,n,jumpI,a] = mvpattern_covcomp(data,z_matrix,'Ac',Ac,'num_iter',2000);
                  
% computing RDM from G matrix
for i=1:nConditions
    for j = 1:nConditions
        if i==j
            r_matrix(i,j) = 1;
        else
            r_matrix(i,j) = G(i,j) / (sqrt(G(i,i)) * sqrt(G(j,j)));
        end
    end
end
end
