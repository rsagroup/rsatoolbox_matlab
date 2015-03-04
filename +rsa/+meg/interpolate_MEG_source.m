function interpolate_MEG_source (Models, userOptions)

targetResolution = userOptions.targetResolution;
smoothingWidth = userOptions.smoothingWidth;
modelNumber = userOptions.modelNumber; 
modelName = spacesToUnderscores(Models(modelNumber).name);
nSubjects = numel(userOptions.subjectNames);
adjacencyMatrix = calculateMeshAdjacency(targetResolution, smoothingWidth, userOptions);

parfor subject = 1:nSubjects
    upsampled_L = [];
    upsampled_R = [];
    thisSubject = userOptions.subjectNames{subject};

    inputPath = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_rMesh_' modelName '_' thisSubject]);
    lh_Vol = mne_read_stc_file1([inputPath,'-lh.stc']); % Pull in data, requires MNE in the search path
    rh_Vol = mne_read_stc_file1([inputPath,'-rh.stc']); % Pull in data, requires MNE in the search path

    [nVertex, nTimePoints] = size(lh_Vol.data);
    
    lh_Vol.data = sortrows([double(lh_Vol.vertices), lh_Vol.data]);
    rh_Vol.data = sortrows([double(rh_Vol.vertices), rh_Vol.data]);
    lh_Vol.data = lh_Vol.data(:,2:end);
    rh_Vol.data = rh_Vol.data(:,2:end);
    
    lh_Vol.data = [lh_Vol.data ; zeros(targetResolution-nVertex,nTimePoints)];
    upsampled_L.vertices = 1:targetResolution;% changed by IZ 03/12 0:targetResolution-1;
    upsampled_L.tstep = lh_Vol.tstep;
    upsampled_L.tmin = lh_Vol.tmin;
    rh_Vol.data = [rh_Vol.data ; zeros(targetResolution-nVertex,nTimePoints)];
    upsampled_R.vertices = 1:targetResolution;% changed by IZ 03/12 0:targetResolution-1;
    upsampled_R.tstep = rh_Vol.tstep;
    upsampled_R.tmin = rh_Vol.tmin;
    
    for vertex = 1:targetResolution
            % Determine which vertexes are within the radius of the currently-picked vertex
            verticesCurrentlyWithinRadius = adjacencyMatrix(vertex,:);
            verticesCurrentlyWithinRadius = verticesCurrentlyWithinRadius(~isnan(verticesCurrentlyWithinRadius));        
            verticesCurrentlyWithinRadius = [vertex, verticesCurrentlyWithinRadius]; 

            for t = 1:nTimePoints
                l = lh_Vol.data(verticesCurrentlyWithinRadius,t);
                r = rh_Vol.data(verticesCurrentlyWithinRadius,t);
                upsampled_L.data(vertex,t) = mean(l(l~=0),1);
                upsampled_R.data(vertex,t) = mean(r(r~=0),1);
            end%for:t
            % Indicate progress every once in a while...
            if mod(vertex, floor(targetResolution/2)) == 0, fprintf('.'); end%if
    end%for:vertex
    
    mne_write_stc_file1([inputPath, '_up-lh.stc'], upsampled_L);   
    mne_write_stc_file1([inputPath, '_up-rh.stc'], upsampled_R); 
end% for:subject
fprintf('Done.\n');