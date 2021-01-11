function distance=fsl_demo()

userOptions=projectOptions_fsl_demo;

[betas,residuals]=rsa.fsl.getDataFromFSL(userOptions,'T05');

partQ = kron(1:length(userOptions.run_names),ones(1,length(userOptions.copes)))';  %runNumber (partition)
partK = repmat(1:length(userOptions.copes),[1 length(userOptions.run_names)])';    %condNumber (CondVec)

u_hat   = rsa.fsl.noiseNormalizeBetaFSL(betas,residuals,partQ);   % Get noise normalised betas
distance  = rsa.distanceLDC(u_hat,partQ,partK);

figure; imagesc(squareform(distance));



end