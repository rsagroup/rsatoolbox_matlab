% MEGSearchlight_source
%
% MEGSearchlight_source(Models, userOptions)
%
% It is based on Su Li's code
%
% Cai Wingfield 3-2010, 9-2010 Su Li updated 3-2012

function [varargout] = MEGSearchlight_source(subjectNumber, Models, userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later
% TODO: This should be passed into the function, and not calculated
% TODO: in-line.  betaCorrespondence.m is something modified by the user,
% TODO: and as such should be called by the Recipe code (also
% TODO: user-modifable), and not by library functions.  This will make it
% TODO: easier in future to separate the library code from a specific
% TODO: recipe format.
tempBetas = userOptions.betaCorrespondence();
subject = userOptions.subjectNames{subjectNumber};
nSubjects = userOptions.nSubjects;

modelNumber = userOptions.modelNumber;
modelName = spacesToUnderscores(Models(modelNumber).name);

if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end

if userOptions.maskingFlag
    MapsFilename = [userOptions.analysisName, '_rMesh_', modelName, '_', subject '_masked'];
else
    MapsFilename = [userOptions.analysisName, '_rMesh_', modelName, '_', subject];
end

promptOptions.functionCaller = 'MEGSearchlight_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-lh.stc']);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-rh.stc']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    
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
        
        % Get masks
        
        % update IZ 03/12
        % assuming every analysis works as a mask, overlaying masks to create a single mask
        nMasks = numel(fieldnames(userOptions.indexMasks));
        indexMasks = userOptions.indexMasks;
        maskNames = fieldnames(userOptions.indexMasks);
        maskIndices=[];
        for mask = 1:nMasks
            thisMaskName = maskNames{mask};
            if strfind(thisMaskName, [lower(chi), 'h'])
                maskIndices = union(maskIndices, indexMasks.(thisMaskName).maskIndices);
                maskIndices = sort(maskIndices(maskIndices <= userOptions.nVertices));
                userOptions.maskIndices.(chi) = maskIndices;
                userOptions.chi = chi;
            end
        end
        timeIndices = userOptions.dataPointsSearchlightLimits;
        
        % Apply masks
        
        if userOptions.nSessions == 1
            maskedMesh = sourceMeshes.(chi)(:, timeIndices(1):timeIndices(2), :); % (vertices, timePointes, conditions)
        else
            maskedMesh = sourceMeshes.(chi)(:, timeIndices(1):timeIndices(2), :, :); % (vertices, timePointes, conditions, sessions)
        end
        [thisSubjectRs.(chi), thisSubjectPs.(chi), searchlightRDMs] = searchlightMapping_MEG_source(maskedMesh, Models, userOptions);
        
        rMetadataStruct = userOptions.STCmetaData;
        pMetadataStruct = userOptions.STCmetaData;
        
        rMetadataStruct.data = thisSubjectRs.(chi)(:,:,modelNumber);
        pMetadataStruct.data = thisSubjectPs.(chi)(:,:,modelNumber);
        
        %% Saving r-maps and p-maps
        outputRFilename = fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_rMesh_' modelName '_' subject ]);
        outputPFilename = fullfile(userOptions.rootPath, 'Maps', modelName,  [userOptions.analysisName '_pMesh_' modelName '_' subject ]);
        if userOptions.maskingFlag
            outputRFilename = [outputRFilename '_masked'];
            outputPFilename = [outputPFilename '_masked'];
        end
        mne_write_stc_file1([outputRFilename '-' lower(chi) 'h.stc'], rMetadataStruct);
        mne_write_stc_file1([outputPFilename '-' lower(chi) 'h.stc'], pMetadataStruct);
        
        %% Saving the searchlight RDMs
        fprintf('Saving data RDMs for combined mask: ');
        filepath = 'searchlightRDMs_';
        if userOptions.maskingFlag
            filepath = [filepath 'masked_'];
        end
        gotoDir(userOptions.rootPath, 'RDMs');
        save('-v7.3', [filepath, userOptions.subjectNames{subjectNumber},'-',lower(chi),'h'], 'searchlightRDMs');
        
        userOptions = rmfield(userOptions, 'maskIndices');
        userOptions = rmfield(userOptions, 'chi');
        clear thisSubjectRs thisSubjectPs pMetadataStruct searchlightRDMs sourceMeshes.(chi) maskedMesh;
        
    end%for:chirality
    
    % Print the elapsed time for this subject
    if chirality == 1
        fprintf('\n\t\t\t\t\t\t\t\t');
    else
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
