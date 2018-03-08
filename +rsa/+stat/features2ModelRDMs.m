function Model=features2ModelRDMs(features,varargin); 
% function RDMs=rsa_features2ModelRDMs({Feature1,'name1',Feature2,'name2'},options); 
% generates a Model structure from features
% for parametric analysis of RDMs 
Opt.dissimilarity='LDC';        % Could be covariance 
import rsa.*; 

Opt=getUserOptions(varargin,Opt,{'dissimilarity'}); 
if (~iscell(features) || ~mod(length(features),2)==0) 
    error('features needs to be a cell array {Feature1,''name1'',Feature2,''name'',...}');
end; 

K = size(features{1},1); 
C = rsa.util.indicatorMatrix('allpairs',[1:K]); 

for i=1:length(features)/2
    Model(i).name=features{i*2}; 
    F   = features{i*2-1}; 
    if (size(F,1)~=K)
        error('all features need to be nConditions x nFeatures'); 
    end; 
    switch (Opt.dissimilarity) 
        case {'LDC','sqEucledian'}; 
            Model(i).IPM = F*F';                        % Previously F*pinv(F'*F)*F': I guess this was to normalise feature vectors?? 
            Model(i).RDM = squareform(diag(C*Model(i).IPM*C')); 
        otherwise 
            error('Don''t know how to get from features to model for dissimilarity other than LDC'); 
    end; 
end; 
