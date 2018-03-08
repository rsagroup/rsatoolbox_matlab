function Model=modelMANOVA(Fa,Fb,varargin);
% function Model=modelMANOVA(Fa,Fb,varargin)
% This function generates the RDM model for a fully-crossed MANOVA design (2-factorial).
% INPUT:
%   Fa: a nConditions column-vector, indicating the assignment of condition k to
%       level A
%   Fb: a nConditions column-vector, indicating the assignment of condition k to
%       level B
% OUTPUT:
%   Model components for the simple main effects, cross-effects, and Interaction
%   The fitting is thought to be applied to LDC distances 
% VARARGIN / USEROPTIONS:
%   Either structure or list of arguments with 
%   'factorNameA', name: Name of the factor A
%   'factorNameB', name: Name of the factor B 
% Joern Diedrichsen 
% joern.diedrichsen@googlemail.com 
% July 2015

Opt.factorNameA='A';               % Factor name A
Opt.factorNameB='B';               % Factor name B
Opt.dissimilarity = 'LDC';         % Dissimilarity measure

% Allow for user options
Opt=rsa.getUserOptions(varargin,Opt,{'factorNameA','factorNameB','dissimilarity'});

% Get sizes of variables, input checking
Ka = length(unique(Fa));
Kb = length(unique(Fb));
K = size(Fa,1);
K_alt = size(Fb,1);
if (K==1 | K_alt==1)
    error('Fa and Fb should be column vectors');
end;
if (K~=K_alt)
    error('Fa and Fb need to have same length');
end;

% make condition matrices
CA=rsa.util.indicatorMatrix('identity',Fa);
CB=rsa.util.indicatorMatrix('identity',Fb);
CI=rsa.util.indicatorMatrix('interaction',[Fa Fb]);

% From this, Calcualte thm into LDC-RDMs
Pairs=rsa.util.indicatorMatrix('allpairs',[1:K]);
Model(1).name      = Opt.factorNameA;
Model(1).IPM       = rsa_vectorizeIPM(CA*CA');
Model(1).RDM       = diag(Pairs*CA*CA'*Pairs')';
Model(1).color     = [1 0 0 ];
        
Model(2).name      = Opt.factorNameB;
Model(2).IPM       = rsa_vectorizeIPM(CB*CB');
Model(2).RDM       = diag(Pairs*CB*CB'*Pairs')';
Model(2).color     = [0 0 1];
        
Model(3).name      = [Opt.factorNameA ' x ' Opt.factorNameB ];
Model(3).IPM       = rsa_vectorizeIPM(CI*CI');
Model(3).RDM       = diag(Pairs*CI*CI'*Pairs')';
Model(3).color     = [1 0 1];

