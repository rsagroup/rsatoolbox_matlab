function [betas,residuals] = getDataFromFSL(userOptions,subjectName)

%load mask
curmaskPath=rsa.util.replaceWildcards(userOptions.maskPath, ...
        '[[subjectName]]',subjectName, ...
        '[[maskName]]',userOptions.maskNames{1});
mask=load_nii(curmaskPath);

%load copes
betas=[];
residuals=[];
for r = userOptions.run_names
    curfeatsPath=rsa.util.replaceWildcards(userOptions.featsPath, ...
        '[[subjectName]]',subjectName, ...
        '[[featPrefix]]',userOptions.featsPrefix, ...
        '[[runName]]',r{1}, ...
        '[[featSuffix]]',userOptions.featsSuffix);
    for c = userOptions.copes
        beta=load_nii(fullfile(curfeatsPath,'stats',['cope' int2str(c{1})' '.nii.gz']));
        betas=cat(2,betas,beta.img(boolean(mask.img)));
    end
    
    res=load_nii(fullfile(curfeatsPath,'stats','res4d.nii.gz'));
    res1=res.img(repmat(boolean(mask.img),[1 1 1 size(res.img,4)]));
    residuals=cat(2,residuals,reshape(res1,[length(res1)/size(res.img,4) size(res.img,4)]));
    %ress=cat(2,ress,res.img(boolean(mask.img)));
end

residuals=residuals';
betas=betas';

end
