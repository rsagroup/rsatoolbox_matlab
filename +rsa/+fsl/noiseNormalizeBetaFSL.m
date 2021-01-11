function [u_hat] = noiseNormalizeBetaFSL(Beta,Res,Partition)

[T,P] = size(Res);
[Q,P] = size(Beta);

Nrun= max(Partition);
partT = kron(1:Nrun,ones(1,T/Nrun))'; 
partQ = Partition;

for i=1:Nrun
    idxT=partT==i;
    idxQ=partQ==i;
    Sw_hat(:,:,i)=rsa.stat.covdiag(Res(idxT,:));       %%% regularize Sw_hat through optimal shrinkage
    u_hat(idxQ,:)=Beta(idxQ,:)*Sw_hat(:,:,i)^(-1/2);   %%% multivariate noise normalization
end;

end