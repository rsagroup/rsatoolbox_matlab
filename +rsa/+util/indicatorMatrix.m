function [X,a]=indicatorMatrix(what,c)
% function X=indicatorMatrix(what,c)
% Makes a indicator Matrix for the classes in c
% INPUT:
%   what   : gives the type of indicator matrix that is needed
%      'identity'       : a regressor for each category 
%      'identity_p'     : a regressor for each category, except for c==0 
%      'reduced'        : GLM-reduced coding (last category is -1 on all indicators)
%      'pairs'          : Codes K category as K-1 pairs
%      'allpairs'       : Codes all K*(K-1) pairs in the same sequence as pdist
%      'allpairs_p'     : Same a allpairs, but ignoring 0 
%      'interaction_reduced' 
%      'hierarchical'   : Codes for a fixed variable underneath a normal factor 
%      'hierarchicalI'  : Codes for a random variable unerneath a normal factor
%   c      : 1xN or Nx1 vector of categories
% OUTPUT:
%   X      : design Matrix
% Joern Diedrichsen 2012

transp=0;
if size(c,1)==1
    c      =  c';
    transp =  1;
end;

[row,col]  = size(c);

for s=1:size(c,2) 
    [a,~,cc(:,s)]=unique(c(:,s)); % Make the class-labels 1-K 
end;

K=max(cc);     % number of classes, assuming numbering from 1...max(c)

switch (what)
    case 'identity'                 % Dummy coding matrix 
        X=zeros(row,K);
        for i=1:K
            X(cc==i,i) = 1;
        end;
    case 'identity_fix'             % Dummy coding matrix with fixed columns 1...K
        X=zeros(row,max(c));
        for i=1:max(c)
            X(c==i,i) = 1;
        end;    
    
    case 'identity_p'               % Dummy coding matrix except for 0: Use this to code simple main effects 
        a(a==0)=[]; 
        K=length(a); 
        X=zeros(row,K);
        for i=1:K
            X(c==a(i),i) = 1;
        end;
        
    case 'reduced'                  % Reduced rank dummy coding matrix (last category all negative) 
        X=zeros(row,K-1);
        for i=1:K-1
            X(cc==i,i) = 1;
        end;
        X(sum(X,2)==0,:)=-1; 
    case 'reduced_p'                % Reduced rank dummy coding matrix, except for 0 
        a(a==0)=[]; 
        K=length(a); 
        X=zeros(row,K-1);
        for i=1:K-1
            X(c==a(i),i) = 1;
        end;
        X(c==a(K),:)=-1; 
        
    case 'pairs'                    % K-1 pairs 
        X=zeros(row,(K-1)); 
        for i=1:K-1
            X(cc==i,i)   =  1./sum(cc==i);
            X(cc==i+1,i) = -1./sum(cc==i+1);
        end;
    case 'allpairs'                 % all possible pairs 
        X=zeros(row,K*(K-1)/2);
        k=1;
        for i=1:K
            for j=i+1:K
                X(cc==i,k) = 1./sum(cc==i);
                X(cc==j,k) = -1./sum(cc==j);
                k         = k+1;
            end;
        end;
    case 'allpairs_p'               % all possible pairs  except for 0 
        a(a==0)=[]; 
        K=length(a); 
        X=zeros(row,K*(K-1)/2);
        k=1;
        for i=1:K
            for j=i+1:K
                X(c==a(i),k) = 1./sum(c==a(i));
                X(c==a(j),k) = -1./sum(c==a(j));
                k         = k+1;
            end;
        end;
    case 'interaction' 
        X1=indicatorMatrix('identity',cc(:,1)); 
        X2=indicatorMatrix('identity',cc(:,2)); 
        for n=1:row
            X(n,:)=kron(X1(n,:),X2(n,:)); 
        end; 
    case 'interaction_reduced' 
        X1=indicatorMatrix('reduced',cc(:,1)); 
        X2=indicatorMatrix('reduced',cc(:,2)); 
        for n=1:row
            X(n,:)=kron(X1(n,:),X2(n,:)); 
        end; 
    case 'hierarchical'             % Allows for a random effects f c(:,2) within each level of c(:,1)
        X1=indicatorMatrix('identity',cc(:,1)); 
        X2=indicatorMatrix('reduced',cc(:,2)); 
        for n=1:row
            X(n,:)=kron(X1(n,:),X2(n,:)); 
        end; 
    case 'hierarchicalI'             % Allows for a random effects f c(:,2) within each level of c(:,1)
        C=size(cc,2); 
        for i=1:C
            X{i}=indicatorMatrix('identity',cc(:,i)); 
        end;
        A=X{end}; 
        for i=C-1:-1:1
            B=[]; 
            for n=1:row
                B(n,:)=kron(X{i}(n,:),A(n,:)); 
            end; 
            A=B; 
        end; 
        X=A; 
    otherwise
        error('no such case');
end;

% Transpose design matrix
if (transp==1)
    X = X';
else
    a=a';
end;
