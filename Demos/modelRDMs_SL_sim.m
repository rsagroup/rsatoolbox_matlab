function Models = modelRDMs_SL_sim()
% 
Models.main_clusters = kron([
			0 0 1 1
			0 0 1 1
			1 1 0 0
			1 1 0 0], ones(10,10));
% Models.random = squareform(pdist(rand(40,40)));

end%function
