%FJ 10/2014
% This function will delte users jobs from the queue 

function flushQ
myCluster = parcluster('CBU_Cluster')
delete(myCluster.Jobs);

end