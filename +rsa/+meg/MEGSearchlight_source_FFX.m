% MEGSearchlight_source
%
% MEGSearchlight_source(Models, userOptions)
%
% It is based on Su Li's code
%
% Cai Wingfield 3-2010, 9-2010 Su Li updated 3-2012

function [varargout] = MEGSearchlight_source_FFX(Models, userOptions)

returnHere = pwd; % We'll come back here later

tempBetas = betaCorrespondence;  

modelNumber = userOptions.modelNumber; 
modelName = spacesToUnderscores(Models(modelNumber).name);
nSubjects = userOptions.nSubjects;

MapsFilename = [userOptions.analysisName, '_rMesh_', modelName, subject];
promptOptions.functionCaller = 'MEGSearchlight_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-lh.stc']);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-rh.stc']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    averageRDMs = [];
    for subjectNumber = 1:nSubjects
        subject = userOptions.subjectNames{subjectNumber};
        
        fprintf('Shining RSA searchlights...\n');
        gotoDir(userOptions.rootPath, 'Maps');
        gotoDir(fullfile(userOptions.rootPath, 'Maps'), modelName); 

        tic;%1
        fprintf(['\tSearching in the source meshes of subject ' num2str(subjectNumber) ' of ' num2str(nSubjects) ':']);

        % Run searchlight on both halves of the brain
        sourceMeshes = MEGDataPreparation_source(subjectNumber, tempBetas, userOptions);

        for chirality = 1:2
            switch chirality
                case 1
                    chi = 'L';
                case 2
                    chi = 'R';
            end%switch:chirality

            % Reduce to temporal RoI
            searchableMeshes.(chi) = sourceMeshes.(chi)(:, userOptions.temporalSearchlightLimits(1):userOptions.temporalSearchlightLimits(2),:,:);

            % Do the searchlight! ZOMG, this takes a while...
            if strcmp(userOptions.groupStats, 'FFX')
                thisSubjectRDMs.(chi) = searchlightMapping_MEG_source_FFX(searchableMeshes.(chi), Models, userOptions); 
                averageRDMs.(chi) = averageRDMs.(chi) + thisSubjectRDMs.(chi);
            else
                [thisSubjectRs.(chi), thisSubjectPs.(chi)] = searchlightMapping_MEG_source_RFX(searchableMeshes.(chi), Models, userOptions); 
                
                % Store results in indexed volumes
                rMeshes.(modelName).(chi) = thisSubjectRs.(chi)(:,:,modelNumber); % r-values for correlation with each model
                pMeshes.(modelName).(chi) = thisSubjectPs.(chi)(:,:,modelNumber);

                clear thisSubjectRs thisSubjectPs;

                % Output filename
                outputRFilename = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_rMesh_' modelName '_' subject '-' lower(chi) 'h.stc']);
                outputPFilename = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_pMesh_' modelName '_' subject '-' lower(chi) 'h.stc']);

                rMetadataStruct = userOptions.STCmetaData;
                pMetadataStruct = userOptions.STCmetaData;

                rMetadataStruct.data = rMeshes.(modelName).(chi);
                pMetadataStruct.data = pMeshes.(modelName).(chi);

        %                 ignore = evalc('gotoDir(userOptions.rootPath, ''Maps'');'); clear ignore;
        %                 gotoDir(fullfile(userOptions.rootPath, 'Maps'), modelName); 

                mne_write_stc_file(outputRFilename, rMetadataStruct);
                mne_write_stc_file(outputPFilename, pMetadataStruct);
                clear  pMetadataStruct rMeshes pMeshes nMeshes;     
            end
            
            clear searchableMeshes;

        end%for:chirality
        
        

        if chirality == 1
            fprintf('\n\t\t\t\t\t\t\t\t');
        else
            % Print the elapsed time for this subject
            t = toc;%1
            fprintf([': [' num2str(ceil(t)) 's]\n']);
        end%if

        %% Save relevant info
    %     timeStamp = datestr(now);
    % 
    %     gotoDir(userOptions.rootPath, 'Details');
    %     cd(fullfile(userOptions.rootPath, 'Details'));
    %     fprintf(['Saving Details to ' fullfile(pwd, DetailsFilename) '\n']);
    %     save(DetailsFilename, 'timeStamp', 'userOptions');
    end
else
    fprintf('Searchlight already applied, skip....\n');
end
% 
% gotoDir(userOptions.rootPath, 'ImageData');
% cd(fullfile(userOptions.rootPath, 'ImageData'));
% fprintf(['Saving Details to ' fullfile(pwd, MetaDataFilename) '\n']);
% save(MetaDataFilename, 'searchlightOptions');

if nargout == 1
    varargout{1} = [];
elseif nargout > 0
    error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started