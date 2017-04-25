function [varargout] = fMRISearchlightModelComparison(models,binaryMasks_nS, betaCorrespondence, userOptions)
%
% fMRISearchlightModelComparison is a function which takes some pre-computed searchlightRDMs,
% and some models and perfoms a searchlight in the data matching to each of the models.
% Saved are native-space r-maps for each model.
%
% [rMaps_sS, maskedSmoothedRMaps_sS] = fMRISearchlightModelComparison(searchlightRDMs,
%                                                 models,
%                                                 userOptions)
%
%        models --- A stack of model RDMs in a structure.
%               models is a [1 nModels] structure with fields:
%                       RDM
%                       name
%                       color
% 		 betaCorrespondence --- The array of beta filenames.
%               betas(condition, session).identifier is a string which referrs
%               to the filename (not including path) of the SPM beta image. (Or,
%               if not using SPM, just something, as it's used to determine the
%               number of conditions and sessions.)
%               Alternatively, this can be the string 'SPM', in which case the
%               SPM metadata will be used to infer this information, provided
%               that userOptions.conditionLabels is set, and the condition
%               labels are the same as those used in SPM.
%
%        userOptions --- The options struct.
%                userOptions.projectName
%               		 A string which is prepended to the saved files.
%		    	   	     This string is specific to the current project
%		  		 userOptions.analysisName
%		         	     A string which is prepended to the saved files.
%		   				 This string is specific to the current analysis.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in fullBrainVols.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Defaults to the fieldnames of the first subject
%                        of binaryMasks_nS.
%                userOptions.structuralsPath
%                        A string which contains the absolute path to the
%                        location of the structural images and the normalisation
%                        warp definition file. It can contain the following
%                        wildcards which would be replaced as indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%
% The following files are saved by this function:
%        userOptions.rootPath/Maps/
%                userOptions.projectName_userOptions.analysisName_fMRISearchlight_Maps.mat
%                        Contains the searchlight statistical maps in struct so
%                        that rMaps_nS.(modelName).(subject).(maskName),
%                        rMaps_sS.(modelName).(subject).(maskName),
%                        maskedSmoothedRMaps_sS.(modelName).(subject).(maskName)
%                        and nMaps_nS.(modelName).(subject).(maskName) contain
%                        the appropriate data.
%        userOptions.rootPath/Details/
%                userOptions.projectName_userOptions.analysisName_fMRISearchlight_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
% this function is dependent on fmriPrepareSearchlightRDMs from which we load:
%
%        searchlightRDMs --- The searchlight RDMs.
% 			    see help rsa.fmri.fMRIPrepareSearchlightRDMs
%
%
% Cai Wingfield 2-2010, 3-2010
% Ian Charest 3-2017 refactored to a 2-step procedure 1 - fmriPrepareSearchlightRDMs 2 - fMRISearchlightModelComparison
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

%% Set defaults and check options struct
if ~isfield(userOptions, 'projectName'), error('fMRISearchlight:NoProjectName', 'ProjectName must be set. See help'); end%if
if ~isfield(userOptions, 'analysisName'), error('fMRISearchlight:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('fMRISearchlight:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'writeOut', 1);
if ~isfield(userOptions, 'voxelSize'), error('fMRISearchlight:NoVoxelSize', 'voxelSize must be set. See help'); end%if


% The analysisName will be used to label the files which are eventually saved.
mapsFilename = sprintf('%s_%s_fMRISearchlight_Maps.mat',userOptions.projectName,userOptions.analysisName);
RDMsFilename = sprintf('%s_fMRISearchlight_RDMs.mat',userOptions.projectName); % the brain RDMs don't change from analysis to analysis.
DetailsFilename = sprintf('%s_%s_fMRISearchlight_Details.mat',userOptions.projectName,userOptions.analysisName);

% load the pre-computed searchlightRDMs, return if no pre-computed searchlightRDMs.
try
    load(fullfile(userOptions.rootPath, 'RDMs',RDMsFilename));
    userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(searchlightRDMs));
    userOptions = setIfUnset(userOptions, 'maskNames', fieldnames(searchlightRDMs.(userOptions.subjectNames{1})));
catch
    error('fMRISearchlightModelComparison:NoSearchlightRDMs , searchlightRDMs not found. Please run fMRIPrepareSearchlightRDMs before model comparison; skipping.');
end

promptOptions.functionCaller = 'fMRISearchlightModelComparison';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', mapsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    
    % Data
    nSubjects = numel(userOptions.subjectNames);
    nMasks = numel(userOptions.maskNames);
    
    searchlightOptions.monitor = false;
    searchlightOptions.fisher = true;
    
    % assing warp flags if writing files
    if userOptions.writeOut==1
        warpFlags.interp = 1;
        warpFlags.wrap = [0 0 0];
        warpFlags.vox = userOptions.voxelSize; % [3 3 3.75]
        warpFlags.bb = [-78 -112 -50; 78 76 85];
        warpFlags.preserve = 0;
    end
    
    fprintf('Shining RSA searchlights...\n');
    
    for subjectNumber = 1:nSubjects % and for each subject...
        
        tic;%1
        
        fprintf(['\t...in the brain of subject ' num2str(subjectNumber) ' of ' num2str(nSubjects) '...\n']);
        
        % Figure out which subject this is
        subject = userOptions.subjectNames{subjectNumber};
        
        
        if userOptions.writeOut==1
            if ischar(betaCorrespondence) && strcmpi(betaCorrespondence, 'SPM')
                betas = getDataFromSPM(userOptions);
            else
                betas = betaCorrespondence;
            end%if:SPM
            
            readFile = replaceWildcards(userOptions.betaPath, '[[subjectName]]', subject, '[[betaIdentifier]]', betas(1,1).identifier);
            subjectMetadataStruct = spm_vol(readFile);
        end
        %		subjectMetadataStruct = spawnSPMStruct;
        
        for maskNumber = 1:nMasks % For each mask...
            
            % Get the mask
            maskName = userOptions.maskNames{maskNumber};
            mask = binaryMasks_nS.(subject).(maskName);
            
            % Do the searchlight! ZOMG, this takes a while...
            [rs, ps] = searchlightModelComparison_fMRI(searchlightRDMs.(subject).(maskName), models, userOptions, searchlightOptions); % ps are from linear correlation p-values, and so aren't too useful here.
            
            for modelNumber = 1:numel(models)
                
                modelName = spacesToUnderscores(models(modelNumber).name);
                
                % Store results in indexed volumes
                rMaps_nS.(modelName).(subject).(maskName) = rs(:,:,:,modelNumber); % r-values for correlation with each model
                
                % Store results in indexed volumes
                pMaps_nS.(modelName).(subject).(maskName) = ps(:,:,:,modelNumber); % p-values for correlation with each model
                
                if userOptions.writeOut==1
                    
                    %% Save native space version
                    
                    % Write the native-space r-map to a file
                    rMapMetadataStruct_nS = subjectMetadataStruct;
                    rMapMetadataStruct_nS.fname = fullfile(userOptions.rootPath, 'Maps', [userOptions.analysisName '_rMap_' maskName '_' modelName '_' subject '.img']);
                    rMapMetadataStruct_nS.descrip =  'R-map';
                    rMapMetadataStruct_nS.dim = size(rMaps_nS.(modelName).(subject).(maskName));
                    
                    gotoDir(userOptions.rootPath, 'Maps');
                    
                    rsa.spm.spm_write_vol(rMapMetadataStruct_nS, rMaps_nS.(modelName).(subject).(maskName));
                    
                    if isfield(userOptions, 'structuralsPath')
                        
                        % Write the native-space mask to a file
                        maskMetadataStruct_nS = subjectMetadataStruct;
                        maskMetadataStruct_nS.fname = fullfile(userOptions.rootPath, 'Maps', [userOptions.analysisName '_nativeSpaceMask_' maskName '_' modelName '_' subject '.img']);
                        maskMetadataStruct_nS.descrip =  'Native space mask';
                        maskMetadataStruct_nS.dim = size(mask);
                        
                        rsa.spm.spm_write_vol(maskMetadataStruct_nS, mask);
                        
                        % Load in common space warp definition
                        % 					wildFiles = replaceWildcards(fullfile(userOptions.structuralsPath, ['*' subject '*_seg_sn.mat']), '[[subjectName]]', subject);
                        wildFiles = replaceWildcards(fullfile(userOptions.structuralsPath, ['*_seg_sn.mat']), '[[subjectName]]', subject);
                        matchingFiles = dir(wildFiles);
                        warpDefFilename = replaceWildcards(fullfile(userOptions.structuralsPath, matchingFiles(1).name), '[[subjectName]]', subject);
                        
                        % Warp and write common space r-maps to disk
                        spm_write_sn(rMapMetadataStruct_nS,warpDefFilename,warpFlags);
                        
                        % Warp and write common space masks to disk
                        spm_write_sn(maskMetadataStruct_nS,warpDefFilename,warpFlags);
                        
                        % Now read them back in
                        
                        % Where are they?
                        [warpedPath_rMap, warpedFile_rMap, warpedExt_rMap, warpedVersion_rMap] = fileparts(rMapMetadataStruct_nS.fname);
                        [warpedPath_mask, warpedFile_mask, warpedExt_mask, warpedVersion_mask] = fileparts(maskMetadataStruct_nS.fname);
                        
                        % Warped versions are prefixed with 'w'
                        warpedFile_rMap = ['w' warpedFile_rMap];
                        warpedFile_mask = ['w' warpedFile_mask];
                        
                        % Read them from the disk
                        rMaps_sS.(modelName).(subject).(maskName) = spm_read_vols(spm_vol(fullfile(warpedPath_rMap, [warpedFile_rMap warpedExt_rMap]))); % sS for standard space
                        mask_sS = spm_read_vols(spm_vol(fullfile(warpedPath_mask, [warpedFile_mask warpedExt_mask])));
                        
                        % Fix the normalisation of the mask
                        maskMetadataStruct_sS = spm_vol(fullfile(warpedPath_rMap, [warpedFile_rMap warpedExt_rMap]));
                        maskMetadataStruct_sS.fname = fullfile(userOptions.rootPath, 'Maps', [userOptions.analysisName '_commonSpaceMask_' maskName '_' modelName '_' subject '.img']);
                        maskMetadataStruct_sS.descrip =  'Common space mask';
                        maskMetadataStruct_sS.dim = size(mask_sS);
                        
                        maskThreshold = 0.01;
                        mask_sS(mask_sS < maskThreshold) = 0;
                        mask_sS(isnan(mask_sS)) = 0;
                        
                        maskMetadataStruct_sS.dim = size(mask_sS);
                        
                        rsa.spm.spm_write_vol(maskMetadataStruct_sS, mask_sS);
                        
                        % Smooth the normalised data
                        
                        % Smoothed versions are prefixed with 's'
                        smoothedWarpedFile_rMap = ['s' warpedFile_rMap];
                        
                        % Smooth it
                        smoothingKernel_fwhm = [10 10 10];
                        spm_smooth(fullfile(warpedPath_rMap, [warpedFile_rMap warpedExt_rMap]), fullfile(warpedPath_rMap, [smoothedWarpedFile_rMap warpedExt_rMap]), smoothingKernel_fwhm);
                        
                        % Read it back in
                        smoothedDataMetadataStruct = spm_vol(fullfile(warpedPath_rMap, [smoothedWarpedFile_rMap warpedExt_rMap]));
                        smoothedData = spm_read_vols(smoothedDataMetadataStruct);
                        
                        % Mask the smoothed data by the sS mask
                        maskedData = smoothedData;
                        maskedData(mask_sS == 0) = NaN;
                        maskedSmoothedRMaps_sS.(modelName).(subject).(maskName) = maskedData;
                        
                        % Write it back to disk
                        maskedDataMetadataStruct_nS = smoothedDataMetadataStruct;
                        maskedDataMetadataStruct_nS.fname = fullfile(userOptions.rootPath, 'Maps', ['msw' userOptions.analysisName '_rMap_' maskName '_' modelName '_' subject '.img']); % 'msw' for 'masked, smoothed, warped'
                        maskedDataMetadataStruct_nS.descrip =  'Masked smoothed normalised data';
                        maskedDataMetadataStruct_nS.dim = size(maskedData);
                        
                        rsa.spm.spm_write_vol(maskedDataMetadataStruct_nS, maskedData);
                        
                    end%if:structuralsPath
                
                end%if:writeOut
                
            end%for:models
            
            clear fullBrainVolumes rs ps ns;
            
            fprintf(':');
            
        end%for:maskNumber
        
        t = toc;%1
        fprintf([' [' num2str(ceil(t)) 's]\n']);
        
    end%for:subjectNumber
    
    %% Save relevant info
    
    timeStamp = datestr(now);
    
    fprintf(['Saving searchlight maps to ' fullfile(userOptions.rootPath, 'Maps', mapsFilename) '\n']);
    gotoDir(userOptions.rootPath, 'Maps');
    if isfield(userOptions, 'structuralsPath')
        save(mapsFilename, 'rMaps_nS', 'rMaps_sS', 'maskedSmoothedRMaps_sS', 'nMaps_nS','pMaps_nS');
    else
        save(mapsFilename, 'rMaps_nS','pMaps_nS','nMaps_nS');
    end%if
    
    fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
    gotoDir(userOptions.rootPath, 'Details');
    save(DetailsFilename, 'timeStamp', 'userOptions');
    
else
    fprintf(['Loading previously saved maps from ' fullfile(userOptions.rootPath, 'Maps', mapsFilename) '...\n']);
    load(fullfile(userOptions.rootPath, 'Maps', mapsFilename));
end%if

if nargout == 2 % native space if 2 variables out
    varargout{1} = rMaps_nS;
    varargout{2} = pMaps_nS;
elseif nargout == 4 % both native space and warped images if 4 variables out
    varargout{1} = rMaps_sS;
    varargout{2} = maskedSmoothedRMaps_sS;
    varargout{3} = rMaps_nS;
    varargout{4} = nMaps_nS;
elseif nargout > 0
    error('0, 2 or 4 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started

end%function


%%%%%%%%%%%%%%%%%%%
%% Sub functions %%
%%%%%%%%%%%%%%%%%%%


function [smm_rs, smm_ps] = searchlightModelComparison_fMRI(searchlightRDMs, models, userOptions, localOptions)

% ARGUMENTS
% searchlightRDMs A nPairWiseComparison x voxel matrix of RDMs.
%
% models		A struct of model RDMs.
%
% userOptions and localOptions
%
% RETURN VALUES
% smm_rs        4D array of 3D maps (x by y by z by model index) of
%               correlations between the searchlight pattern similarity
%               matrix and each of the model similarity matrices.
%
% smm_ps        4D array of 3D maps (x by y by z by model index) of p
%               values computed for each corresponding entry of smm_rs.
%
% Based on Niko Kriegeskorte's searchlightMapping_RDMs.m
%
% Additions by Cai Wingfield 2-2010:
% 	- Now skips points in the searchlight where there's only one voxel inside.
% 	- Now takes a userOptions struct for the input parameters.
%
% Refactored to a 2-step model by Ian Charest 3-2017

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%% Get parameters
monitor = localOptions.monitor;

% Prepare models
modelRDMs_ltv = permute(unwrapRDMs(vectorizeRDMs(models)), [3 2 1]);

% Other data
volSize_vox=size(searchlightRDMs);
nPairWiseComparisons = volSize_vox(1);
volSize_vox = volSize_vox(2:end);
nModelRDMs=size(modelRDMs_ltv,1);

%% similarity-graph-map the volume with the searchlight
smm_bestModel=nan(volSize_vox);
smm_ps=nan([volSize_vox,nModelRDMs]);
smm_rs=nan([volSize_vox,nModelRDMs]);

%% find the indices of all nan RDMs

nVox_mappingMask_request = prod(volSize_vox);
tsearchlightRDMs = reshape(searchlightRDMs,[nPairWiseComparisons,nVox_mappingMask_request]);

nonNanIs = find(not(all(isnan(tsearchlightRDMs)))); clear searchlightRDMs

if monitor
    h_progressMonitor=progressMonitor(1, nVox_mappingMask_request,  'Similarity-graph-mapping...');
end

%% THE BIG LOOP! %%

for cMappingVoxI=nonNanIs
    
    if mod(cMappingVoxI,1000)==0
        if monitor
            progressMonitor(cMappingVoxI, nVox_mappingMask_request, 'Searchlight model comparison ...', h_progressMonitor);
            %                 cMappingVoxI/nVox_mappingMask_request
        else
            fprintf('.');
        end%if
    end%if
    
    [x,y,z]=ind2sub(volSize_vox,cMappingVoxI);
    
    % extract the searchlight RDM
    
    thisRDM = tsearchlightRDMs(:,cMappingVoxI)';
    
    if not(all(isnan(thisRDM)))
    
        try
            [rs, ps] = corr(thisRDM', modelRDMs_ltv', 'type', 'Spearman', 'rows', 'pairwise');
        catch
            [rs, ps] = corr(thisRDM', modelRDMs_ltv, 'type', 'Spearman', 'rows', 'pairwise');
        end%try

        if localOptions.fisher
            for i = 1:numel(rs)
                rs(i) = fisherTransform(rs(i));
            end%for:i
        end%if

        %	[ignore, bestModelI] = max(rs);

        %    smm_bestModel(x,y,z) = bestModelI;
        smm_ps(x,y,z,:) = ps;
        smm_rs(x,y,z,:) = rs;
    end
end%for:cMappingVoxI

%% END OF THE BIG LOOP! %%
fprintf('\n');
if monitor
    close(h_progressMonitor);
end

%% visualize
if monitor
    aprox_p_uncorr=0.001;
    singleModel_p_crit=aprox_p_uncorr/nModelRDMs; % conservative assumption model proximities nonoverlapping
    
    smm_min_p=min(smm_ps,[],4);
    smm_significant=smm_min_p<singleModel_p_crit;
    
    vol=map2vol(mask);
    vol2=map2vol(mask);
    
    colors=[1 0 0
        0 1 0
        0 1 1
        1 1 0
        1 0 1];
    
    for modelRDMI=1:nModelRDMs
        vol=addBinaryMapToVol(vol, smm_significant&(smm_bestModel==modelRDMI), colors(modelRDMI,:));
        % 		vol2=addBinaryMapToVol(vol2, smm_bestModel==modelRDMI, colors(modelRDMI,:));
    end
    
    showVol(vol);
    
    % 	vol2 = vol2*0.1;
    % 	vol2(vol<1) = vol(vol<1);
    %
    % 	showVol(vol2);
    
end%if

end%function
