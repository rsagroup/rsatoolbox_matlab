% This function computes group statistics across sliding time windows using
% Random Effects Analysis.
% This will print and save an uncorrected t/r-map and a corrected t/r-map
% thresholded based on cluster statistics. Each cluster will have a mass
% and a corresponding p-value.
% Input:    userOptions, Models

% Based on scripts by Li Su
% Written by IZ 03/13 updated FJ 03/14



function RFX_slidingTimeWindow(userOptions, Models)

close all;
returnHere = pwd; % We'll come back here later
modelNumber = userOptions.modelNumber;
modelName = Models(modelNumber).name;
if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end


output_path = fullfile(userOptions.rootPath, 'Results', 'RandomEffects');
promptOptions.functionCaller = 'RFX_slidingTimeWindow';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(output_path, [modelName '-' userOptions.maskNames{numel(userOptions.maskNames)} '-r.xls']);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
    rdms_path =fullfile(userOptions.rootPath,'RDMs',[userOptions.analysisName '_' modelName '_dataRDMs_sliding_time_window']);
    
    if ~exist(output_path,'dir')
        mkdir(output_path);
    end
    
    fprintf('Loading all data RDMs... ');
    try
        load(rdms_path);
        disp('Done!');
    catch
        error('Cannot load data RDMs file.')
    end
    
    nMasks = size(allRDMs,1);
    nTimePoints = size(allRDMs,2);
    nSubjects = userOptions.nSubjects;
    
    disp('Performing random effects permutation test...')
    
    for mask=1:nMasks
        thisMask = userOptions.maskNames{mask};
        
        disp('Testing (model, RoI) RDM pairs for significance of similarity...');
        
        %% computing correlation across subjects
         parfor timeWindow = 1:nTimePoints
            
            RDMs = squeeze(allRDMs(mask,timeWindow,:))'; % changing format of RDMs to mask, subject,session
            RDMs = averageRDMs_subjectSession(RDMs, 'session');
            
            modelRDM = Models(modelNumber).RDM;
            modelRDM_vec = vectorizeRDM(modelRDM);
            if userOptions.partial_correlation
                control_for_modelRDMs = [];
                for m = 1:size(userOptions.partial_modelNumber,2)
                    control_for_modelRDMs = [control_for_modelRDMs;vectorizeRDM(Models(userOptions.partial_modelNumber{m}).RDM)];
                end
            end
            for subject=1:nSubjects
                if userOptions.partial_correlation
                    [r(subject,1,timeWindow) p(subject,timeWindow)] = partialcorr(vectorizeRDM(RDMs(1,subject).RDM)',modelRDM_vec',control_for_modelRDMs','type',userOptions.distanceMeasure,'rows','pairwise');
                else
                    [r(subject,1,timeWindow) p(subject,timeWindow)] = corr(vectorizeRDM(RDMs(1,subject).RDM)',modelRDM_vec','type',userOptions.distanceMeasure,'rows','pairwise');
                end
            end
        end
        
        %% computing cluster based t/r statistic
        fpmin = 1;
        fpmax = 1;
        tpmin = 1;
        tpmax = size(r,3);
        perm_num = userOptions.significanceTestPermutations;
        test = 1; % test=1 for a 1 sample t-test, test=2 for a paired t-test
        pval = userOptions.primaryThreshold;
        data1 = r;
        data2 = zeros(size(r));
        tmapFlag = userOptions.tmap;
        [clust_stats_pos, clust_stats_neg, base_map, null_distribution] = ...
            permutation_cluster_test_2dtfr_func(data1, ...
            data2,fpmin,fpmax,tpmin,tpmax,perm_num,test,pval,tmapFlag);
        
        if userOptions.tmap
            threshold = tinv(pval,(userOptions.nSubjects-1));
            thresh_map = ((base_map>= abs(threshold)) .* base_map) + ((base_map <= threshold) .* base_map);
            what_map = 't';
        else
            thresh_map = base_map;
            thresh_map(median(p)' > pval) = 0;
            threshold = min(base_map(median(p)' < pval))- 0.0001;
            if isempty(threshold)
                disp('Warning: there is no threshold set');
                threshold=0;
            end
            what_map = 'r';
        end
        
        all_clusters_pos{mask} = clust_stats_pos;
        %all_clusters_neg{mask} = clust_stats_neg;
        
        percent = 0.05;
        null_distribution = sort(null_distribution);
        cluster_level_threshold = null_distribution(ceil(size(null_distribution,2)*(1-percent)));

        %causing error in display                
%         disp('Plotting null distribution...');
%         figure(mask);
%         hist(null_distribution,100);
%         h = findobj(gca, 'Type', 'patch');
%         set(h, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0 0 0]);
%         hold on;
%         yLimits = get(gca, 'YLim');
%         plot(repmat(cluster_level_threshold,1,yLimits(1,2)), 1:yLimits(1,2), ':','color', 'red', 'LineWidth', 2 );
%         text(double(cluster_level_threshold), yLimits(1,2)/2, [' \leftarrow cluster level threshold: ' num2str(cluster_level_threshold,4)], 'FontSize', 10);   
%         title([thisMask ': Null distribution across ' num2str(userOptions.significanceTestPermutations) ' permutations (rfx)']);
%         saveas(figure(mask),fullfile(userOptions.rootPath, 'Results', 'RandomEffects', [thisMask '_null-distribution']),'fig');
%  
%         
        %% saving files
        fprintf(['Saving uncorrected ' what_map '-map... ']);
        xlswrite(fullfile(output_path, [modelName '-' thisMask '-uncorrected_' what_map '.xls']), base_map);
        disp('Done!');
        fprintf(['Saving corrected ' what_map '-map... ']);
        xlswrite(fullfile(output_path, [modelName '-' thisMask '-corrected_' what_map '.xls']), thresh_map);
        disp('Done!');
        fprintf('Saving r values for all subjects...')
        xlswrite(fullfile(output_path, [modelName '-' thisMask '-r' '.xls']), squeeze(r));
        disp('Done!');
        fprintf('Saving p values for all subjects...')
        xlswrite(fullfile(output_path, [modelName '-' thisMask '-p' '.xls']), squeeze(p));
        disp('Done!');
        fprintf('Saving null distribution...')
        xlswrite(fullfile(userOptions.rootPath, 'ImageData', [modelName '-' thisMask '-rfx-nulldistribution.xls']), null_distribution);
        disp('Done!');
        
        if ~exist(fullfile(output_path, 'ClusterStats'),'dir')
            mkdir(fullfile(output_path, 'ClusterStats'));
        end
        xlswrite(fullfile(output_path, 'ClusterStats', [modelName '-' thisMask '-cluster_stats-pos-' what_map '.xls']), clust_stats_pos);
        % xlswrite([modelName '-' thisMask '-cluster_stats-neg'], clust_stats_neg);
        disp('Done!');
        
    end
else
    fprintf('Permutation already performed, skip....\n');
end

