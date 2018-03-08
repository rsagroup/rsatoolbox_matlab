%  modelRDMs is a user-editable function which specifies the models which
%  brain-region RDMs should be compared to, and which specifies which kinds of
%  analysis should be performed.
%
%  Models should be stored in the "Models" struct as a single field labeled
%  with the model's name (use underscores in stead of spaces).
%  
%  Cai Wingfield 11-2009

function Models = modelRDMs_demo2()

Models.allSeparate = kron([
			0 1 1 1
			1 0 1 1
			1 1 0 1
			1 1 1 0], ones(16,16));

Models.mainClusters = kron([
			0 0 1 1
			0 0 1 1
			1 1 0 0
			1 1 0 0], ones(16,16));
        
Models.leftCluster = kron([
			0 2 1 1
			2 0 1 1
			1 1 0 0
			1 1 0 0], ones(16,16));
        
Models.rightCluster = kron([
			0 0 1 1
			0 0 1 1
			1 1 0 2
			1 1 2 0], ones(16,16));
        

% Models.N = kron([
% 			0 1 1 1
% 			1 0 1 1
% 			1 1 0 1
% 			1 1 1 0], ones(16,16));
        
% Models.bad_prototype = [kron([0 .5; .5 0], ones(16,16)) ones(32,32); ones(32,32), 1-eye(32,32)];
Models.random = squareform(pdist(rand(64,64)));

% Models.noStructure = ones(64,64);
% Models.noStructure(logical(eye(64)))=0;

end%function
