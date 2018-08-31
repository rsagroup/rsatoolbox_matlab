function n2f=surfing_nodeidxs2faceidxs(f)
% finds the mapping from vertices to face indices.
% 
% N2F=NODEIDXS2FACEIDXS(F) 
% INPUT:
%   F:     Px3 node indices for P faces (base1)
% OUTPUT: 
%   N2F:   QxR mapping matrix for Q nodes (Q==MAX(F(:))), if R is the
%          maximum number of faces that contain a single node.
%
% Property of N2F: 
%   For certain J: F(I,J)==N  <==>  for certain K: N2F(N,K)==I.
%
% NNO May 2010, updated Sep 2010
%
% See also SURFING_INVERTMAPPING

if size(f,1) ~= 3
    error('expected 3xP face matrix');
end

n2f=surfing_invertmapping(f');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD CODE 
%(NNO: replaced with the above as of Sep 2010)
% 
% 
% nv=max(f(:)); % number of vertices
% mp=zeros(nv,18); %map of vertices to faces; allocate plenty of space
% p=0;
% for k=1:3 % each row in the face matrix
%     im=surfing_invertmapping(f(k,:)'); % see which nodes are contained in the faces k-th column
%     w=size(im,2); % maximum number of duplicates in this column
%     mp(:,p+(1:w))=im; %add where there is still space
%     p=p+w;
% end
% 
% fperrow=sum(mp>0,2); %faces per row
% maxfperrow=max(fperrow); %max number of faces per row (typically 6)
% 
% n2f=zeros(nv,maxfperrow); % a smaller map with as few zeros as possible
% freepos=ones(nv,1); % first free position in each row
% 
% for ci=1:size(mp,2) % every column in mp
%     mpi=mp(:,ci); % element in column ci
%     msk=mpi>0; % nonzero elements in this column
%     for di=1:size(mp,2) % see where we can put the nonzero elements
%         m=freepos==di & msk; % free position and non-zero?
%         if sum(m)>0 
%             n2f(m,di)=mpi(m); % store elements
%             freepos(m)=freepos(m)+1; %increment free position for these elements
%             msk=msk & ~m; % ensure not to set values again
%         end
%     end
% end

% [just a testing function to ensure that this function behaves)
% fis=floor(rand(1,1000)*nf)+1;
% for fi=fis % face index
%     for j=1:3
%         vi=f(fi,j); % vertex idx
%         fs=mp2(vi,:);
%         if sum(fi==fs)==0
%             error('not found in faces %d: %d', fi, vi);
%         end
%     end
% end
% fprintf('OK!!\n');
        
    
    