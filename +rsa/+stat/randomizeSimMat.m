function randomizedSimMat=randomizeSimMat(simMat)

[v,h]=size(simMat);

if v~=h
    % convert from upper triangular form
    simMat=squareform(simMat);
end

scrambledIs=randperm(size(simMat,1));

randomizedSimMat=simMat(scrambledIs,scrambledIs);

if v~=h
    % convert to upper triangular form
    randomizedSimMat=squareform(randomizedSimMat);
end


