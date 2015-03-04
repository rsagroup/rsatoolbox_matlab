% This function permutes between subjects to simulate a null distribution of
% maximum stats.
%
% Li Su 10-02-2012

function RFX_permutation(Models, userOptions)

returnHere = pwd; % We'll come back here later

nSubjects = userOptions.nSubjects;
modelNumber = userOptions.modelNumber;
modelName = spacesToUnderscores(Models(modelNumber).name);

if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end

MapsFilename = ['perm-', num2str(userOptions.significanceTestPermutations), '_', modelName, '_t_map']; %Li Su edited 11/2012

promptOptions.functionCaller = 'MEGFindCluster_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-lh.stc']);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Maps', modelName, [MapsFilename, '-rh.stc']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    fprintf('Permuting (RFX) ...');
%     if ~userOptions.maskingFlag
        for subjectNumber = 1:nSubjects
            subject = userOptions.subjectNames{subjectNumber};
            inputFilename = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_rMesh_' modelName '_' subject]);
            if userOptions.maskingFlag
                inputFilename = [inputFilename , '_masked'];
            end
            MEGDataStcL = mne_read_stc_file1([inputFilename, '-lh.stc']);
            MEGDataStcR = mne_read_stc_file1([inputFilename, '-rh.stc']);
            MEGDataVolL = single(MEGDataStcL.data);
            MEGDataVolR = single(MEGDataStcR.data);
            all_rho(subjectNumber,:,:) = [MEGDataVolL; MEGDataVolR];
        end
        
        numberOfPermutation = userOptions.significanceTestPermutations;
        [numberOfVertex, numberOfTimePoints] = size(MEGDataStcL.data);
        
        [h,p,ci,stats] = ttest(all_rho);
        t_value = squeeze(stats.tstat);
        t_lh = t_value(1:numberOfVertex,:);
        t_rh = t_value(numberOfVertex+1:numberOfVertex * 2,:);
        lh_Vol = MEGDataStcL;
        rh_Vol = MEGDataStcR;
        
        t_lh(isnan(t_lh)) = 0;
        t_rh(isnan(t_rh)) = 0;
        lh_Vol.data = t_lh;
        rh_Vol.data = t_rh;
        
        
        outputFilename = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_tMesh_' modelName '_allSubjects']);
        if userOptions.maskingFlag
           outputFilename = [outputFilename , '_masked'];
        end
        mne_write_stc_file1([outputFilename, '-lh.stc'], lh_Vol);
        mne_write_stc_file1([outputFilename, '-rh.stc'], rh_Vol);
        observed_Vol_L = lh_Vol;
        observed_Vol_R = rh_Vol;
        
        max_t_value = zeros(1,numberOfPermutation);
        
        parfor perm = 1:numberOfPermutation
            rho = all_rho;
            simulated_Vol_L = lh_Vol;
            simulated_Vol_R = rh_Vol;
            
            if mod(perm, floor(numberOfPermutation/40)) == 0, fprintf('\b.'); end%if
            
            for subjectNumber = 1:nSubjects
                %toss a coin
                toss = (rand(1) > 0.5)*2-1;
                rho(subjectNumber,:,:) = squeeze(all_rho(subjectNumber,:,:)).*toss;
            end
            
            [h,p,ci,stats] = ttest(rho);
            t_value = squeeze(stats.tstat);
            
            t_lh = t_value(1:numberOfVertex,:);
            t_rh = t_value(numberOfVertex+1:numberOfVertex * 2,:);
            simulated_Vol_L.data = t_lh;
            simulated_Vol_R.data = t_rh;
            
            
            outputFilename = fullfile(userOptions.rootPath, 'Maps', modelName, ['perm-' num2str(perm) '_' modelName '_t_map']);
            mne_write_stc_file1([outputFilename, '-lh.stc'], simulated_Vol_L);
            mne_write_stc_file1([outputFilename, '-rh.stc'], simulated_Vol_R);
            
            max_t_value(perm) = max(max(t_value));
        end
        
        percent = 0.05; % Update from userOptions.primaryThreshold; IZ 03,12
        t_distribution = sort(max_t_value);
        
        vertex_level_threshold = t_distribution(ceil(size(t_distribution,2)*(1-percent)));
        
        fprintf('\n Writng results corrected at whole brain level using permutation but without using clustering method...\n');
        
        gotoDir(userOptions.rootPath, 'Results');
        outputFileName_sig = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName '_significant_vertex']);
        if userOptions.maskingFlag
              outputFileName_sig = [outputFileName_sig , '_masked'];
        end
        observed_Vol_L.data(observed_Vol_L.data<vertex_level_threshold) = 0;
        observed_Vol_R.data(observed_Vol_R.data<vertex_level_threshold) = 0;
        
        mne_write_stc_file1([outputFileName_sig, '-lh.stc'], observed_Vol_L);
        mne_write_stc_file1([outputFileName_sig, '-rh.stc'], observed_Vol_R);
    
else
    fprintf('Already done permutation, Skip...');
end

cd(returnHere); % And go back to where you started
