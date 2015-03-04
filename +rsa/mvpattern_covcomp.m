function [G,h,u,l,n,jumpI,a]=mvpattern_covcomp(y,Z,varargin)
%mvpattern_covcomp: estimate random-effects variance component coefficients
% Usage: [G,h,u,l,n,jumpI,a]=mvpattern_covcomp(y,Z,varargin);
% Estimates the variance coefficients of the model described in:
% Diedrichsen, Ridgway, Friston & Wiestler (2011).
% y_n = Z u_n + e,
%         u ~ (b, G)
%                 G = A A'
%
% y: N x P observations
% Z: N x Q random effects matrix
%
% N: numbers of observations (trials)
% P: numbers of experimental units (voxels)
% Q: number of random effects
%
% VARARGIN:
%   'num_iter'     : Maximal number of iterations
%   'Ac'           : Cell array (Hx1) of components of factor loading 
%                    matrix A = sum(h_m A_m), used to form G = A A'
%   'h0'           : Starting values for the parameters (Hx1 vector)
%                    (otherwise uses starting guess based on Laird & Lange)
%   'TolL'         : Tolerance of the likelihood (l-l'), where l' is
%                    projected likelihood 
%   'accel_method' : Acceleration method used
%                    'none'   : unaccelerated EM algorithm 
%                    'Aitken' : Aitken acceleration with projected jump
%                               based on observed convergence 
%   'meanS'        : Remove the mean for each pattern component (a)
%                    (Logical flag, true by default)
%
% OUTPUT:
%   G     : variance-covariance matrix
%   h     : coefficients
%   u     : hidden patterns components 
%   l     : likelihood
%   n     : number of iterations until convergence
%   jumpI : record of jumps in convergence 
%   a     : Fixed effect means for the pattern components
% 
% Examples:
% See mva_component_examples
%  
% See also: mva_component_examples, spm_reml, spm_reml_sc, spm_sp_reml
% Where spm_* are from the SPM software, http://www.fil.ion.ucl.ac.uk/spm
%
% Copyright 2011 Joern Diedrichsen, j.diedrichsen@ucl.ac.uk

[N,P]=size(y);
[N2,Q]=size(Z);
if N2 ~= N
    error('Mismatched numbers of rows in data (%d) and design (%d)', N, N2)
end

num_iter=600;       % Maximal number of iterations
Ac={};
TolL=0.00001;            % Tolerance on Likelihood
accel_method='Aitken';
meanS=1;                    % Mean subtract
h0=[];
% Variable argument otions
vararginoptions(varargin, ...
    {'num_iter','Ac','TolL','accel_method','meanS','h0'});

% Intialize the Model structure
H=length(Ac);
Cc = cell(H, 1);
xC=zeros(Q*Q,H);
for i=1:H
    Cc{i}=Z*Ac{i}; % Calculate regression coefficients for full path C
    xC(:,i)=Ac{i}(:);
end;

% Precompute cross-terms for M-step
CcCc = cell(H, H);
% for i=1:H
%     for j=i:H
%         CcCc{i,j}=Cc{j}'*Cc{i}; % Calculate cross terms
%     end;
% end;

% If necessary, subtract the fixed effects estimates (a)
if (meanS)
    a=pinv(Z'*Z*P)*Z'*sum(y,2);
    r=y-repmat(Z*a,1,P);
    rr=r*r';                        % This is Suffient stats 1 (S1)
else
    b=zeros(size(Z,2),1);
    r=y-repmat(Z*b,1,P);
    rr=r*r';                        % S1
end;

% preallocate arrays for speed
l=zeros(1,num_iter);
h=zeros(H+1,num_iter); % h(H+1) = Sigma2
delta_h=zeros(H+1,num_iter);

% Provide initial guess from Laird, Lange and Dempster
if (isempty(h0))
    u=pinv(Z'*Z)*Z'*r;
    h(H+1,1)=(sum(sum(y.*y,2))-a'*Z'*sum(y,2)-sum(sum((u'*Z').*r')))/(P*N);
    D=u*u'/P-h(H+1,1)*pinv(Z'*Z)/P;
    hD=diag(diag(D).^0.5);
    % Starting values for constrained estimates
    h(1:H,1)=(xC'*xC)\(xC'*hD(:));
else 
    h(1:H,1)=h0; 
end; 

% Initialize
n=1;
jump=1;
jumpI=nan(1, num_iter);
jumpI(1)=1;
diffL=inf;

% Iterate
while (n<num_iter && diffL>TolL)
    
    % Build up current C-matrix
    C=zeros(N,Q);
    for i=1:H
        C=C+Cc{i}*h(i,n);
    end;
    
    % Estep
    V=h(H+1,n)*eye(N)+C*C';
    Wr=V\r;
    WC=V\C;
    v=C'*Wr;
    rv=r*v';                            % <rv'> suffienct stats (S2)
    vv=(v*v')+P*eye(Q)-P*C'*WC;          % <vv'> sufficient statistics (S3)
    l(n)=-P/2*(log(det(V)))-0.5*sum(sum(Wr.*r,2));
    
    % Check if likelihood decreased on the last iteration
    if (n>1 && l(n)<l(n-1) || h(H+1,n)<0)
        
        % Check if last iteration was a jump
        if (~jumpI(n)==1) % It wasn't: that's bad it should not decrease
            warning('mvpattern_covcomp:EMdecrease', ...
                'EM decreased by %g', l(n-1)-l(n));
            diffL=0;
            if (l(n-1)-l(n)>0.001)
                % If it is only a small decrease, it may be rouding error
                disp('likelihood decreased!');
            end;
        else
            % Last step was a jump: So go undo the jump and redo the E-step
            if n>1, n=n-1; end
            
            % Build up current C-matrix
            C=zeros(N,Q);
            for i=1:H
                C=C+Cc{i}*h(i,n);
            end;
            
            % Estep
            V=h(H+1,n)*eye(N)+C*C';
            Wr=V\r;
            WC=V\C;
            v=C'*Wr;
            rv=r*v';
            vv=(v*v')+P*eye(Q)-P*C'*WC; % <vv'> sufficient statistics
            l(n)=-P/2*(log(det(V)))-0.5*sum(sum(Wr.*r,2));
        end;
    end;
    
    % Mstep: Constrained regression
    % This part is the biggest time sink and would likely profit highly
    % from being rewritten as a mex file. It uses the following trace-trick
    % trace(A*B)=sum(sum(A.*B')), which is faster to compute.
    COV = zeros(H,1);
    VA = zeros(H, H);
    for i = 1:H
        COV(i,1)=sum(sum(Cc{i}.*rv));           % trace-trick
%         VA(i,i)=sum(sum(CcCc{i,i}.*vv));        % trace-trick
          VA(i,i)=sum(sum((Cc{i}'*Cc{i}).*vv));
        for j = i+1:H                           % off-diagonal terms
%             VA(i,j)=sum(sum(CcCc{i,j}.*vv));    % trace-trick
            VA(i,j)=sum(sum((Cc{j}'*Cc{i}).*vv));
            VA(j,i)=VA(i,j);                    % asumme symmetry
        end;
    end;
    h(1:H,n+1) = VA\COV;
    
    
    % Based on the new h-parameters, build up C
    C=zeros(N,Q);
    for i=1:H
        C=C+Cc{i}*h(i,n+1);
    end;
    
    % Compute the residuals and estimate sigma_e
    R=rr-C*rv';
    h(H+1,n+1)=1/(N*P)*trace(R);
    
    % Track the change in parameters
    delta_h(:,n+1)=h(:,n+1)-h(:,n);
    
    % Track the change in likelihood as a stopping criterion.
    % Do not abort for small steps in liklihood, but only when the
    % real likelihood is close to the estimated maximal likelihood.
    % This prevents stopping of the iteration when there is only slow
    % progress (see McLaughlan & Krishnan, 1997. The EM alogithm and
    % Extensions)
    if (n-jump>2)                           % Check convergence by
        Rdl=(l(n)-l(n-1))/(l(n-1)-l(n-2));  % Ratio of differences in l
        lA=l(n-1)+1./(1-Rdl)*(l(n)-l(n-1)); % Predicted maximal likelihood
        diffL=lA-l(n);                      % Estimated deviation
    end;
    
    % prepare next iteration
    n=n+1;
    jumpI(n)=0;
    
    % See if a jump can be done based on recent progress
    if (strcmp(accel_method,'Aitken'))
        if (n-jump>3)
            lambda_hat=mean(delta_h(:,n)./delta_h(:,n-1));
            if (lambda_hat<1)
                h(:,n+1)=h(:,n-1)+1./(1-lambda_hat)*(h(:,n)-h(:,n-1));
                l(n)=l(n-1);
                n=n+1;
                jumpI(n)=1;
                jump=n;
            end
        end
    end
end
jumpI(n+1:end)=[]; % trim down from num_iter length if diffL converged

% Now build the G-matrix of random effects
% First build the full A matrix, then build G
% This was a problem before, because I ignored the crossterms
% sum(Ac_i*Ac_i'*h_i*h_i) != sum(Ac_i*h_i)*sum(Ac_i*h_i)'
h=h(:,1:n-1);
l=l(1:n-1);
A=zeros(Q,Q);
for i=1:H
    A=A+Ac{i}*h(i,end);
end;
G=A*A';
u=G*Z'*Wr;

function vararginoptions(options,allowed_vars,allowed_flags)
% function vararginoptions(options,allowed_vars,allowed_flags);
% Deals with variable argument in   
% INPUTS
%   options: cell array of a argument list passed to a function
%   allowed_vars: Variables that can be set 
%   allowed_flags: Flags that can be set 
%  vararginoptions assigns the value of the option to a variable with the
%  name option (in called workspace).
%  Flags are set to one:
% EXAMPLE:
%   the option-string 'var1',4,'var2',10,'flag'
%   causes the var1 and var2 to be set to 4 and 10 and flag to 1
%   if allowedvars are not given, all variables are allowed
% Joern Diedrichsen 
% v1.0 9/13/05
checkflags=1;
checkvars=1;
if nargin<2
    checkvars=0;
end;
if nargin<3
    checkflags=0;
end;

c=1;
while c<=length(options)
    a=[];
    if ~ischar(options{c})
        error('Options must be strings on argument %d',c);
    end;
    if checkflags
        a=find(strcmp(options{c},allowed_flags), 1);
    end;
    if ~isempty(a)
        assignin('caller',options{c},1);
        c=c+1;
    else
        if checkvars
            a=find(strcmp(options{c},allowed_vars), 1);
            if (isempty(a))
                error(['unknown option:' options{c}]);
            end;
        end;
        if (c==length(options))
            error('Option %s must be followed by a argument',options{c});
        end;
        assignin('caller',options{c},options{c+1});
        c=c+2;
    end;
end;
