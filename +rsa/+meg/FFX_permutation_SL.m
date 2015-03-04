% This function permutes between subjects to simulate a null distribution of
% maximum stats.

% Based on RFX script by Su Li
% Isma Zulfiqar 11/2012 Updated IZ 04/13

function FFX_permutation (Models, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later
nSubjects = numel(userOptions.subjectNames);
modelNumber = userOptions.modelNumber;
modelName = spacesToUnderscores(Models(modelNumber).name);

MapsFilename = ['perm-', userOptions.significanceTestPermutations, '_', modelName, '_r_map'];

promptOptions.functionCaller = 'FFX_permutation';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-lh.stc']);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-rh.stc']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    disp('Permuting (FFX) ...');

    %%  computing average subject RDMs
    for chirality = 1:2
        switch chirality
            case 1
                chi = 'L';
            case 2
                chi = 'R';
        end % switch: chilarity
        
        % computing combined mask
        nMasks = numel(fieldnames(userOptions.indexMasks));
        indexMasks = userOptions.indexMasks;
        masks = fieldnames(userOptions.indexMasks);
        maskIndices.(chi)=[];
        for mask = 1:nMasks
            thisMask = masks{mask};
            if strfind(thisMask,[lower(chi),'h'])
                maskIndices.(chi) = union(maskIndices.(chi),indexMasks.(thisMask).maskIndices);
                maskIndices.(chi) = sort(maskIndices.(chi)(maskIndices.(chi) <= userOptions.nVertices));
            end
        end
        nVertices = length(maskIndices.(chi));
        
        fprintf(['Computing Average Subject RDMs for ' (chi) ' hemisphere...']);
        
        
        MapsFilename = 'averaged_searchlightRDMs_masked_';

        promptOptions.functionCaller = 'FFX_permutation';
        promptOptions.defaultResponse = 'S';
        promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'RDMs', [MapsFilename, 'lh.mat']);
        promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'RDMs', [MapsFilename, 'rh.mat']);

        overwriteFlag = overwritePrompt(userOptions, promptOptions);
        
        if overwriteFlag
        % computes nanmean
            for subjectNumber = 1:nSubjects
                filepath = 'searchlightRDMs_';
                if userOptions.maskingFlag
                    filepath = [filepath 'masked_'];
                end
                subjectRDMsFile = fullfile(userOptions.rootPath, 'RDMs', [filepath 'subject_', userOptions.subjectNames{subjectNumber} '-' lower(chi) 'h']);
                subjectRDMs = load(subjectRDMsFile);

                for v=1:nVertices % vertices
                    nTimePoints = length(fieldnames(subjectRDMs.searchlightRDMs.(['v_' num2str(maskIndices.(chi)(v))])));
                    for t=1:nTimePoints % time points
                        averageSubjectRDMs.(chi).(['v_' num2str(maskIndices.(chi)(v))]).([...
                            't_' num2str(t)]).RDM(subjectNumber,:) = vectorizeRDM...
                            (subjectRDMs.searchlightRDMs.(['v_' num2str(maskIndices.(chi)(v))]).([...
                            't_' num2str(t)]).RDM);

                        if subjectNumber == nSubjects
                            averageSubjectRDMs.(chi).(['v_' num2str(maskIndices.(chi)(v))]).([...
                                't_' num2str(t)]).RDM = nanmean(averageSubjectRDMs.(chi).([...
                                'v_' num2str(maskIndices.(chi)(v))]).(['t_' num2str(t)]).RDM,1);
                        end
                    end % timepoints
                end % vertices
                fprintf('.');
            end % subjects
            save('-v7.3',fullfile(userOptions.rootPath,'RDMs',['averaged_' filepath  lower(chi) 'h']), 'averageSubjectRDMs');
            clear subjectRDMs;
            disp(' Done!');
        else 
            filepath = 'searchlightRDMs_';
            if userOptions.maskingFlag
                filepath = [filepath 'masked_'];
            end
            load(promptOptions.checkFiles(chirality).address);
            nTimePoints = length(fieldnames(averageSubjectRDMs.(chi).(['v_' num2str(maskIndices.(chi)(1))])));
        end
        
        % preparing models
        modelNumber = userOptions.modelNumber;
        modelRDMs_utv = vectorizeRDMs(Models(1,modelNumber));
        modelRDMs_utv = modelRDMs_utv.RDM;
                
        % observed correlation
        temp = zeros(userOptions.targetResolution, nTimePoints);
        fprintf('Computing observed correlation...');
        for i = 1:nVertices
            vertex = maskIndices.(chi)(i);
            rdms = averageSubjectRDMs.(chi).(['v_' num2str(vertex)]);
            parfor t = 1:nTimePoints
                rdm = rdms.(['t_' num2str(t)]).RDM;
                r = corr(vectorizeRDM(rdm)', modelRDMs_utv', 'type', ...
                    userOptions.distanceMeasure, 'rows', 'pairwise');
                temp(vertex, t) = r;
            end
        end
        observed_r = squeeze(temp);
        disp(' Done!');
        
        % saving file to stc format
        observed_Vol.(chi) = userOptions.STCmetaData;
        observed_Vol.(chi).data = observed_r;
        
        outputFilename = fullfile(userOptions.rootPath, 'Maps', modelName, ...
            [userOptions.analysisName '_rMesh_' modelName '_allSubjects']);
        if userOptions.maskingFlag
            outputFilename = [outputFilename '_masked'];
        end
        mkdir(fullfile(userOptions.rootPath, 'Maps'),modelName);
        mne_write_stc_file1([outputFilename, '-', lower(chi), 'h.stc'], observed_Vol.(chi));
        
        clear averageSubjectRDMs;% rather than storing both simultaneously in ram, load them independently
        clear temp rdms rdm observed_r; 
    end % chirality
    
    %% permuting for generating different models and computing correlation
    
    numberOfPermutation = userOptions.significanceTestPermutations;
    max_r_value = zeros(1,numberOfPermutation);
    
    simulated_Vol_lh = userOptions.STCmetaData;
    simulated_Vol_lh.data = zeros(userOptions.targetResolution, nTimePoints);
    simulated_Vol_rh = userOptions.STCmetaData;
    simulated_Vol_rh.data = zeros(userOptions.targetResolution, nTimePoints);
    
    fprintf('Permuting...');
    for perm = 1:numberOfPermutation
        if mod(perm, floor(numberOfPermutation/20)) == 0, fprintf('\b.'); end%if
        %disp(['perm' num2str(perm)]);
        
        % preparing models
        modelRDMs = randomizeSimMat(Models(1,modelNumber).RDM);
        modelRDMs_utv = squeeze(unwrapRDMs(vectorizeRDMs(modelRDMs)));
        
        offset = 0;
        
        for chirality = 1:2
            switch chirality
                case 1
                    chi = 'L';
                case 2
                    chi = 'R';
            end % switch: chilarity
            
            % compute corr with averageSubjectsRDM
            load(fullfile(userOptions.rootPath,'RDMs',['averaged_' filepath  lower(chi) 'h']));
            
            simulated_r = zeros(userOptions.targetResolution*2, nTimePoints);
            for i = 1:length(maskIndices.(chi))
                vertex = maskIndices.(chi)(i);
                temp = zeros(1,nTimePoints);
                rdms = averageSubjectRDMs.(chi).(['v_' num2str(vertex)]);
                parfor t = 1:nTimePoints
                    rdm = rdms.(['t_' num2str(t)]).RDM;
                    r = corr(vectorizeRDM(rdm)', modelRDMs_utv', 'type', ...
                        userOptions.distanceMeasure, 'rows', 'pairwise');
                    temp(t) = r;
                end % for t
                simulated_r(vertex + offset,:) = temp;
                
            end % for vertex
           offset = userOptions.targetResolution; 
           
           clear averageSubjectRDMs temp rdms;
           fprintf('.');
        end % chirality
        
        outputFilename = fullfile(userOptions.rootPath, 'Maps', modelName, ['perm-' num2str(perm) '_' modelName '_r_map']);
        
        simulated_Vol_lh.data = simulated_r(1:userOptions.targetResolution,:);
        mne_write_stc_file1([outputFilename,'-lh.stc'], simulated_Vol_lh);
        offset = userOptions.targetResolution;
        
        simulated_Vol_rh.data = simulated_r(userOptions.targetResolution + 1:userOptions.targetResolution*2,:);
        mne_write_stc_file1([outputFilename,'-rh.stc'], simulated_Vol_rh);
   
        max_r_value(perm) = max(max(simulated_r));
        
        clear simulated_r;
        
    end % for perm
    disp(' Done!');
    
    percent = 0.05; % Update IZ 03/13
    r_distribution = sort(max_r_value);
    
    vertex_level_threshold = r_distribution(ceil(size(r_distribution,2)*(1-percent)));
    
    disp('Writng results corrected for both hemispheres using permutation but without using clustering method...');
    
    gotoDir(userOptions.rootPath, 'Results');
    outputFileName_sig = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName '_significant_vertex_ffx']);
    if userOptions.maskingFlag
        outputFileName_sig = [outputFileName_sig '_masked'];
    end
    observed_Vol.L.data(observed_Vol.L.data<vertex_level_threshold) = 0;
    observed_Vol.R.data(observed_Vol.R.data<vertex_level_threshold) = 0;
    
    mne_write_stc_file1([outputFileName_sig, '-lh.stc'], observed_Vol.L);
    mne_write_stc_file1([outputFileName_sig, '-rh.stc'], observed_Vol.R);
    
else
    fprintf('Already done permutation, Skip...');
end
%matlabpool close
cd(returnHere); % And go back to where you started

%%  old code
% %%
%     nVertices = size(averageSubjectRDMs.(chi), 3);
%     nTimepoints = size(averageSubjectRDMs.(chi), 4);
%
%     averagePermutationRs = nan(nVertices, nTimepoints, nModels);
%     averagePermutationPs = nan(nVertices, nTimepoints, nModels);
%
%     for vertex = 1:nVertices
%         for t = 1:nTimepoints
%             for model = 1:nModels
%                 [averagePermutationRs(vertex, t, model) averagePermutationPs(vertex, t, model)] ...
%                     = testRDMrelatedness_vectorisedRandomization( ...
%                         averageSubjectRDMs.(chi)(:,:,vertex,t), ...
%                         Models(model).RDM, ...
%                         userOptions.distanceMeasure, ...
%                         userOptions.significanceTestPermutations ...
%                     );
%             end%for:model
%         end%for:t
%         if mod(vertex, 500) == 0
%             fprintf('.');
%         end%if
%     end%for:vertex
%
%     %% Save the p and r meshes in MNE format
%
%     % Read sample metadata structure
%     readPath = replaceWildcards(userOptions.betaPath, '[[betaIdentifier]]', tempBetas(1, 1).identifier, '[[subjectName]]', userOptions.subjectNames{1}, '[[LR]]', lower(chi));
%     metadataStruct = mne_read_stc_file(readPath);
%
%     % Replace the metadata
%     %metadataStruct.tmin = userOptions.temporalSearchlightLimits(1) / 1000; % in seconds
%     metadataStruct.tstep = metadataStruct.tstep * userOptions.temporalSearchlightResolution;
%
%     % Now copy for rs and for ps for each model
%     rMetadataStruct = metadataStruct;
%     pMetadataStruct = metadataStruct;
%
%     for model = 1:nModels
%
%         modelName = spacesToUnderscores(Models(model).name);
%
%         % Get the right data
%         rMetadataStruct.data = squeeze(averagePermutationRs(:,:,model));
%         pMetadataStruct.data = squeeze(averagePermutationPs(:,:,model));
%
%         % Get the right filenames
%         rFilename = [userOptions.analysisName '_sujectAveraged_rMesh_' modelName '-' lower(chi) 'h.stc'];
%         pFilename = [userOptions.analysisName '_sujectAveraged_pMesh_' modelName '-' lower(chi) 'h.stc'];
%
%         % Write the file
%         cd(fullfile(userOptions.rootPath, 'Maps'));
%         mne_write_stc_file(rFilename, rMetadataStruct);
%         mne_write_stc_file(pFilename, pMetadataStruct);
%
%     end%for:model
%
% end%for:chirality
%
% t = toc;%2
% fprintf(['[' num2str(ceil(t)) 's]\n']);
% end
