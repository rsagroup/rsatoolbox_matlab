% This function compute the null distribution from the results of the clustering
% algorithm. It creates brain maps (.stc) of significant (and marginal significant)
% clusters over time. It also generates reports in .txt files, which show
% the information of significant clusters and their p-values. Finally, a
% histogram of the null distribution is ploted with vertex and cluster
% level thresholds displayed on the sreen.
%
% created by Li Su, last update 11-2012

function get_cluster_p_value(Models,userOptions)

modelNumber = userOptions.modelNumber;
modelName = spacesToUnderscores(Models(modelNumber).name);

if userOptions.partial_correlation
    modelName = [modelName, '_partialCorr'];
end

%percent = userOptions.primaryThreshold;
tStep = userOptions.temporalSearchlightResolution;

% nullDistribution_lh = [];
% nullDistribution_rh = [];

% if not(userOptions.maskingFlag)
nullDistribution = [];
thisTimeWindow = userOptions.temporalSearchlightLimits;
if exist(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution.csv']),'file')
    nullDistribution = csvread(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution.csv']));
else
    for i = 1:userOptions.jobSize:userOptions.significanceTestPermutations
        csvFileName = fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution' num2str(i) '.csv']);
        nullDistribution = [nullDistribution,csvread(csvFileName)];
        if not(userOptions.debug), delete(csvFileName); end
    end
    nullDistribution = reshape(nullDistribution, 1, numel(nullDistribution));
    csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution.csv']), nullDistribution);
end

% nullDistribution_lh = reshape(nullDistribution_lh, 1, numel(nullDistribution_lh));
% nullDistribution_rh = reshape(nullDistribution_rh, 1, numel(nullDistribution_rh));

observedFileName = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_tMesh_' modelName '_allSubjects']);
if userOptions.maskingFlag
    observedFileName = [observedFileName, '_masked'];
end
lh_Vol = mne_read_stc_file1([observedFileName,'-lh.stc']); % Pull in data, requires MNE in the search path
rh_Vol = mne_read_stc_file1([observedFileName,'-rh.stc']); % Pull in data, requires MNE in the search path
lh_r_vec = reshape(lh_Vol.data,1,size(lh_Vol.data,1)*size(lh_Vol.data,2));
rh_r_vec = reshape(rh_Vol.data,1,size(rh_Vol.data,1)*size(rh_Vol.data,2));

lh_cluster = mne_read_stc_file1([observedFileName,'_cluster-lh.stc']); % Pull in data, requires MNE in the search path
rh_cluster = mne_read_stc_file1([observedFileName,'_cluster-rh.stc']); % Pull in data, requires MNE in the search path
lh_cluster_vec = reshape(lh_cluster.data,1,size(lh_cluster.data,1)*size(lh_cluster.data,2));
rh_cluster_vec = reshape(rh_cluster.data,1,size(rh_cluster.data,1)*size(rh_cluster.data,2));

%% compute the cluster p values
if isnan(userOptions.primaryThreshold)
    % find the vertex level threshold of t value by sorting the t
    % distribution and select the cut off point of x%
    percent = 0.05;
    r_distribution = sort([lh_r_vec,rh_r_vec]);
    vertex_level_threshold = r_distribution(ceil(size(r_distribution,2)*(1-percent)));
else
    pval = userOptions.primaryThreshold;
    vertex_level_threshold = abs(tinv(pval,(userOptions.nSubjects-1)));
end


cluster_stats = 0;
if cluster_stats == 0 % using cluster mass as the maximum stats
    numberOfCluster_left = max(lh_cluster_vec);
    lh_cluster_size = [0,1];
    rh_cluster_size = [0,1];
    for i = 1:numberOfCluster_left
        lh_cluster_size(i,1) = sum(lh_r_vec(lh_cluster_vec == i)-vertex_level_threshold);
        lh_cluster_size(i,2) = length(find(nullDistribution > lh_cluster_size(i,1)))/size(nullDistribution,2);
    end
    numberOfCluster_right = max(rh_cluster_vec);
    for i = 1:numberOfCluster_right
        rh_cluster_size(i,1) = sum(rh_r_vec(rh_cluster_vec == i)-vertex_level_threshold);
        rh_cluster_size(i,2) = length(find(nullDistribution > rh_cluster_size(i,1)))/size(nullDistribution,2);
    end
end
if cluster_stats == 1
    load cluster_stats;
    nullDistribution_space = struct2cell(permutation);
    nullDistribution_space = cell2mat(squeeze(nullDistribution_space(1,:,:)));
    nullDistribution_time = struct2cell(permutation);
    nullDistribution_time = cell2mat(squeeze(nullDistribution_time(3,:,:)));
    
    number_of_cluster_lh = max(max(lh_cluster.data));
    for j = 1:number_of_cluster_lh
        [cluster_row, cluster_col] = find(lh_cluster.data == j);
        if not(isempty(cluster_col))
            cluster_duration(j) = max(cluster_col) - min(cluster_col) + 1;
            for t = min(cluster_col):max(cluster_col)
                temp_data = lh(:,t);
                temp_cluster = lh_cluster.data(:,t);
                clusterMassByTime(t) = sum(temp_data(temp_cluster == j)-vertex_level_threshold);
            end
            cluster_mass_by_time(j)  = mean(clusterMassByTime);
        else
            cluster_duration(j) = 0;
            cluster_mass_by_time(j) = 0;
        end
        p_space = length(find(nullDistribution_space > cluster_mass_by_time(j)))/size(nullDistribution_space,1);
        p_time = length(find(nullDistribution_time > cluster_duration(j)))/size(nullDistribution_time,1);
        lh_cluster_size(j,1) = (p_space+p_time)/2;
        lh_cluster_size(j,2) = min(p_space,p_time);
    end
    
    number_of_cluster_rh = max(max(rh_cluster.data));
    for j = 1:number_of_cluster_rh
        [cluster_row, cluster_col] = find(rh_cluster.data == j);
        if not(isempty(cluster_col))
            cluster_duration(j) = max(cluster_col) - min(cluster_col) + 1;
            for t = min(cluster_col):max(cluster_col)
                temp_data = rh(:,t);
                temp_cluster = rh_cluster.data(:,t);
                clusterMassByTime(t) = sum(temp_data(temp_cluster == j)-vertex_level_threshold);
            end
            cluster_mass_by_time(j)  = mean(clusterMassByTime);
        else
            cluster_duration(j) = 0;
            cluster_mass_by_time(j) = 0;
        end
        
        p_space = length(find(nullDistribution_space > cluster_mass_by_time(j)))/size(nullDistribution_space,1);
        p_time = length(find(nullDistribution_time > cluster_duration(j)))/size(nullDistribution_time,1);
        rh_cluster_size(j,1) = (p_space+p_time)/2;
        rh_cluster_size(j,2) = min(p_space,p_time);
    end
    
end

csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_left-p-values.csv']),lh_cluster_size);
csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_right-p-values.csv']),rh_cluster_size);

%% display significant clusters
significant_cluster_lh = find(lh_cluster_size(:,2) <= 0.05);
significant_cluster_rh = find(rh_cluster_size(:,2) <= 0.05);
marginally_significant_cluster_lh = find(lh_cluster_size(:,2) <= 0.1);
marginally_significant_cluster_rh = find(rh_cluster_size(:,2) <= 0.1);

marginally_significant_cluster_lh = marginally_significant_cluster_lh(~ismember(marginally_significant_cluster_lh,significant_cluster_lh));
marginally_significant_cluster_rh = marginally_significant_cluster_rh(~ismember(marginally_significant_cluster_rh,significant_cluster_rh));

outputFileName_sig = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName '_significant_cluster']);
outputFileName_mar_sig = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName '_mar_significant_cluster']);
if userOptions.maskingFlag
    outputFileName_sig = [outputFileName_sig '_masked'];
    outputFileName_mar_sig = [outputFileName_mar_sig '_masked'];
end
temp_lh = lh_cluster;
temp_rh = rh_cluster;

temp_lh.data(~ismember(temp_lh.data,marginally_significant_cluster_lh)) = 0;
temp_lh.data(ismember(temp_lh.data,marginally_significant_cluster_lh)) = 1;
temp_rh.data(~ismember(temp_rh.data,marginally_significant_cluster_rh)) = 0;
temp_rh.data(ismember(temp_rh.data,marginally_significant_cluster_rh)) = 1;

temp_vol_lh = lh_Vol;
temp_vol_rh = rh_Vol;

temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;

temp_vol_lh.tstep = tStep ./ 1000;
temp_vol_rh.tstep = tStep ./ 1000;

mne_write_stc_file1([outputFileName_mar_sig, '-lh.stc'], temp_vol_lh);
mne_write_stc_file1([outputFileName_mar_sig, '-rh.stc'], temp_vol_rh);

temp_lh = lh_cluster;
temp_rh = rh_cluster;

temp_lh.data(~ismember(temp_lh.data,significant_cluster_lh)) = 0;
temp_lh.data(ismember(temp_lh.data,significant_cluster_lh)) = 1;
temp_rh.data(~ismember(temp_rh.data,significant_cluster_rh)) = 0;
temp_rh.data(ismember(temp_rh.data,significant_cluster_rh)) = 1;

temp_vol_lh = lh_Vol;
temp_vol_rh = rh_Vol;

temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;

temp_vol_lh.tstep = tStep ./ 1000;
temp_vol_rh.tstep = tStep ./ 1000;

mne_write_stc_file1([outputFileName_sig, '-lh.stc'], temp_vol_lh);
mne_write_stc_file1([outputFileName_sig, '-rh.stc'], temp_vol_rh);

%% Generate report of the cluster level stats
report =cell(1,8);
for i = 1:size(significant_cluster_lh,1)
    report(i,1) = {significant_cluster_lh(i)}; % cluster number
    report(i,2) = {(find(sum(lh_cluster.data ==significant_cluster_lh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
    report(i,3) = {(find(sum(lh_cluster.data ==significant_cluster_lh(i))>0, 1, 'last' )-1)*tStep+thisTimeWindow(1)}; % end time in ms
    report(i,4) = {(find(sum(lh_cluster.data ==significant_cluster_lh(i))==max(sum(lh_cluster.data ==significant_cluster_lh(i))),1,'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
    temp_lh = lh_cluster;
    temp_lh.data(~ismember(temp_lh.data,significant_cluster_lh(i))) = 0;
    temp_lh.data(ismember(temp_lh.data,significant_cluster_lh(i))) = 1;
    temp_vol_lh = lh_Vol;
    temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
    [x,y]=find(temp_vol_lh.data==max(max(temp_vol_lh.data)));
    report(i,5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
    report(i,6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
    report(i,7) = {lh_cluster_size(significant_cluster_lh(i),1)}; % cluster stats, i.e. cluster mass
    report(i,8) = {lh_cluster_size(significant_cluster_lh(i),2)}; % cluster level p value after the correction
end
for i = 1:size(marginally_significant_cluster_lh,1)
    report(i+size(significant_cluster_lh,1),1) = {marginally_significant_cluster_lh(i)}; % cluster number
    report(i+size(significant_cluster_lh,1),2) = {(find(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
    report(i+size(significant_cluster_lh,1),3) = {(find(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))>0, 1, 'last' )-1)*tStep+thisTimeWindow(1)}; % end time in ms
    report(i+size(significant_cluster_lh,1),4) = {(find(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))==max(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))),1, 'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
    temp_lh = lh_cluster;
    temp_lh.data(~ismember(temp_lh.data,marginally_significant_cluster_lh(i))) = 0;
    temp_lh.data(ismember(temp_lh.data,marginally_significant_cluster_lh(i))) = 1;
    temp_vol_lh = lh_Vol;
    temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
    [x,y]=find(temp_vol_lh.data==max(max(temp_vol_lh.data)));
    report(i+size(significant_cluster_lh,1),5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
    report(i+size(significant_cluster_lh,1),6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
    report(i+size(significant_cluster_lh,1),7) = {lh_cluster_size(marginally_significant_cluster_lh(i),1)}; % cluster stats, i.e. cluster mass
    report(i+size(significant_cluster_lh,1),8) = {lh_cluster_size(marginally_significant_cluster_lh(i),2)}; % cluster level p value after the correction
end

% alternate file writing in table format
fileID = fopen(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_left-report.txt']),'w');
if fileID~=-1 && ~isempty(report)
    fprintf(fileID, '%5i %30s %36s %36s %44s %55s %50s %50s %55s\n','', 'c1: cluster number','c2: start time in ms', 'c3: end time in ms', 'c4: peak time with maximum space', 'c5: peak vertex number with maximum fit, i.e. r value',...
        'c6:peak time with maximum fit, i.e. r value', 'c7: cluster stats, i.e. cluster mass', 'c8: cluster level p value after the correction');
    totalClusters = size(report);
    for i = 1:totalClusters(1)
        fprintf(fileID, '%5i %25g %36g %36g %44g %55g %50g %50g %55g\n', i, flipud(cell2mat(report(i,:))));
    end
end

% if ~isempty(report)
%     dlmwrite(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_left-report.txt']),cell2mat(report),'\t');
% end

report =cell(1,8);
for i = 1:size(significant_cluster_rh,1)
    report(i,1) = {significant_cluster_rh(i)}; % cluster number
    report(i,2) = {(find(sum(rh_cluster.data ==significant_cluster_rh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
    report(i,3) = {(find(sum(rh_cluster.data ==significant_cluster_rh(i))>0, 1, 'last')-1)*tStep+thisTimeWindow(1)+thisTimeWindow(1)}; % end time in ms
    report(i,4) = {(find(sum(rh_cluster.data ==significant_cluster_rh(i))==max(sum(rh_cluster.data ==significant_cluster_rh(i))),1, 'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
    temp_rh = rh_cluster;
    temp_rh.data(~ismember(temp_rh.data,significant_cluster_rh(i))) = 0;
    temp_rh.data(ismember(temp_rh.data,significant_cluster_rh(i))) = 1;
    temp_vol_rh = rh_Vol;
    temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;
    [x,y]=find(temp_vol_rh.data==max(max(temp_vol_rh.data)));
    report(i,5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
    report(i,6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
    report(i,7) = {rh_cluster_size(significant_cluster_rh(i),1)}; % cluster stats, i.e. cluster mass
    report(i,8) = {rh_cluster_size(significant_cluster_rh(i),2)}; % cluster level p value after the correction
end
for i = 1:size(marginally_significant_cluster_rh,1)
    report(i+size(significant_cluster_rh,1),1) = {marginally_significant_cluster_rh(i)}; % cluster number
    report(i+size(significant_cluster_rh,1),2) = {(find(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
    report(i+size(significant_cluster_rh,1),3) = {(find(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))>0, 1, 'last')-1)*tStep+thisTimeWindow(1)}; % end time in ms
    report(i+size(significant_cluster_rh,1),4) = {(find(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))==max(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))),1, 'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
    temp_rh = rh_cluster;
    temp_rh.data(~ismember(temp_rh.data,marginally_significant_cluster_rh(i))) = 0;
    temp_rh.data(ismember(temp_rh.data,marginally_significant_cluster_rh(i))) = 1;
    temp_vol_rh = rh_Vol;
    temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;
    [x,y]=find(temp_vol_rh.data==max(max(temp_vol_rh.data)));
    report(i+size(significant_cluster_rh,1),5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
    report(i+size(significant_cluster_rh,1),6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
    report(i+size(significant_cluster_rh,1),7) = {rh_cluster_size(marginally_significant_cluster_rh(i),1)}; % cluster stats, i.e. cluster mass
    report(i+size(significant_cluster_rh,1),8) = {rh_cluster_size(marginally_significant_cluster_rh(i),2)}; % cluster level p value after the correction
end

% alternate file writing in table format
fileID = fopen(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_right-report.txt']),'w');
if fileID~=-1 && ~isempty(report)
    fprintf(fileID, '%5i %30s %36s %36s %44s %55s %50s %50s %55s\n','', 'c1: cluster number','c2: start time in ms', 'c3: end time in ms', 'c4: peak time with maximum space', 'c5: peak vertex number with maximum fit, i.e. r value',...
        'c6:peak time with maximum fit, i.e. r value', 'c7: cluster stats, i.e. cluster mass', 'c8: cluster level p value after the correction');
    totalClusters = size(report);
    for i = 1:totalClusters(1)
        fprintf(fileID, '%5i %25g %36g %36g %44g %55g %50g %50g %55g\n', i, flipud(cell2mat(report(i,:))));
    end
end

% if ~isempty(report)
%     dlmwrite(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_right-report.txt']),cell2mat(report),'\t');
% end

fprintf(['vertex level threshold = ', num2str(vertex_level_threshold), '\n']);
nullDistribution = sort(nullDistribution);
cluster_level_threshold = nullDistribution(ceil(size(nullDistribution,2)*(1-0.05)));
fprintf(['cluster level threshold = ', num2str(cluster_level_threshold), '\n']);
% hist(nullDistribution,0:max(nullDistribution)/100:max(nullDistribution));

% else
%     nMasks = numel(userOptions.maskNames);
%     %for mask = 1:2:nMasks % crude way (assumption: every mask is paired lh and rh)
%         nullDistribution = [];
%         %thisMask_lh = dashToUnderscores(userOptions.maskNames{mask});
%         %[thisMask hemisphere] = strtok(thisMask_lh,'_');
%         thisTimeWindow = userOptions.dataPointsSearchlightLimits;
%
%         for i = 1:userOptions.jobSize:userOptions.significanceTestPermutations
%             csvFileName = fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution' num2str(i) '.csv']);
%             nullDistribution = [nullDistribution,csvread(csvFileName)];
%             if not(userOptions.debug), delete(csvFileName); end
%         end
%
%         nullDistribution = reshape(nullDistribution, 1, numel(nullDistribution));
%         csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_null_distribution.csv']), nullDistribution);
%         % nullDistribution_lh = reshape(nullDistribution_lh, 1, numel(nullDistribution_lh));
%         % nullDistribution_rh = reshape(nullDistribution_rh, 1, numel(nullDistribution_rh));
%
%         observedFileName = fullfile(userOptions.rootPath, 'Maps', modelName, [userOptions.analysisName '_tMesh_' modelName '_allSubjects']);
%         lh_Vol = mne_read_stc_file1([observedFileName,'-lh.stc']); % Pull in data, requires MNE in the search path
%         rh_Vol = mne_read_stc_file1([observedFileName,'-rh.stc']); % Pull in data, requires MNE in the search path
%         lh_r_vec = reshape(lh_Vol.data,1,size(lh_Vol.data,1)*size(lh_Vol.data,2));
%         rh_r_vec = reshape(rh_Vol.data,1,size(rh_Vol.data,1)*size(rh_Vol.data,2));
%
%         lh_cluster = mne_read_stc_file1([observedFileName,'_cluster-lh.stc']); % Pull in data, requires MNE in the search path
%         rh_cluster = mne_read_stc_file1([observedFileName,'_cluster-rh.stc']); % Pull in data, requires MNE in the search path
%         lh_cluster_vec = reshape(lh_cluster.data,1,size(lh_cluster.data,1)*size(lh_cluster.data,2));
%         rh_cluster_vec = reshape(rh_cluster.data,1,size(rh_cluster.data,1)*size(rh_cluster.data,2));
%
%         %% compute the cluster p values
%         if isnan(userOptions.primaryThreshold)
%             % find the vertex level threshold of t value by sorting the t
%             % distribution and select the cut off point of x%
%             percent = 0.05;
%             r_distribution = sort([lh_r_vec,rh_r_vec]);
%             vertex_level_threshold = r_distribution(ceil(size(r_distribution,2)*(1-percent)));
%         else
%             pval = userOptions.primaryThreshold;
%             vertex_level_threshold = abs(tinv(pval,(userOptions.nSubjects-1)));
%         end
%
%
%         cluster_stats = 0;
%         if cluster_stats == 0 % using cluster mass as the maximum stats
%             numberOfCluster_left = max(lh_cluster_vec);
%             lh_cluster_size = [0,1];
%             rh_cluster_size = [0,1];
%             for i = 1:numberOfCluster_left
%                 lh_cluster_size(i,1) = sum(lh_r_vec(lh_cluster_vec == i)-vertex_level_threshold);
%                 lh_cluster_size(i,2) = length(find(nullDistribution > lh_cluster_size(i,1)))/size(nullDistribution,2);
%             end
%             numberOfCluster_right = max(rh_cluster_vec);
%             for i = 1:numberOfCluster_right
%                 rh_cluster_size(i,1) = sum(rh_r_vec(rh_cluster_vec == i)-vertex_level_threshold);
%                 rh_cluster_size(i,2) = length(find(nullDistribution > rh_cluster_size(i,1)))/size(nullDistribution,2);
%             end
%         end
%         if cluster_stats == 1
%             load cluster_stats;
%             nullDistribution_space = struct2cell(permutation);
%             nullDistribution_space = cell2mat(squeeze(nullDistribution_space(1,:,:)));
%             nullDistribution_time = struct2cell(permutation);
%             nullDistribution_time = cell2mat(squeeze(nullDistribution_time(3,:,:)));
%
%             number_of_cluster_lh = max(max(lh_cluster.data));
%             for j = 1:number_of_cluster_lh
%                 [cluster_row cluster_col] = find(lh_cluster.data == j);
%                 if not(isempty(cluster_col))
%                     cluster_duration(j) = max(cluster_col) - min(cluster_col) + 1;
%                     for t = min(cluster_col):max(cluster_col)
%                         temp_data = lh(:,t);
%                         temp_cluster = lh_cluster.data(:,t);
%                         clusterMassByTime(t) = sum(temp_data(temp_cluster == j)-vertex_level_threshold);
%                     end
%                     cluster_mass_by_time(j)  = mean(clusterMassByTime);
%                 else
%                     cluster_duration(j) = 0;
%                     cluster_mass_by_time(j) = 0;
%                 end
%                 p_space = length(find(nullDistribution_space > cluster_mass_by_time(j)))/size(nullDistribution_space,1);
%                 p_time = length(find(nullDistribution_time > cluster_duration(j)))/size(nullDistribution_time,1);
%                 lh_cluster_size(j,1) = (p_space+p_time)/2;
%                 lh_cluster_size(j,2) = min(p_space,p_time);
%             end
%
%             number_of_cluster_rh = max(max(rh_cluster.data));
%             for j = 1:number_of_cluster_rh
%                 [cluster_row cluster_col] = find(rh_cluster.data == j);
%                 if not(isempty(cluster_col))
%                     cluster_duration(j) = max(cluster_col) - min(cluster_col) + 1;
%                     for t = min(cluster_col):max(cluster_col)
%                         temp_data = rh(:,t);
%                         temp_cluster = rh_cluster.data(:,t);
%                         clusterMassByTime(t) = sum(temp_data(temp_cluster == j)-vertex_level_threshold);
%                     end
%                     cluster_mass_by_time(j)  = mean(clusterMassByTime);
%                 else
%                     cluster_duration(j) = 0;
%                     cluster_mass_by_time(j) = 0;
%                 end
%
%                 p_space = length(find(nullDistribution_space > cluster_mass_by_time(j)))/size(nullDistribution_space,1);
%                 p_time = length(find(nullDistribution_time > cluster_duration(j)))/size(nullDistribution_time,1);
%                 rh_cluster_size(j,1) = (p_space+p_time)/2;
%                 rh_cluster_size(j,2) = min(p_space,p_time);
%             end
%
%         end
%
%         csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_left-p-values.csv']),lh_cluster_size);
%         csvwrite(fullfile(userOptions.rootPath, 'ImageData', [userOptions.analysisName '_' modelName '_right-p-values.csv']),rh_cluster_size);
%
%         %% display significant clusters
%         significant_cluster_lh = find(lh_cluster_size(:,2) <= 0.05);
%         significant_cluster_rh = find(rh_cluster_size(:,2) <= 0.05);
%         marginally_significant_cluster_lh = find(lh_cluster_size(:,2) <= 0.1);
%         marginally_significant_cluster_rh = find(rh_cluster_size(:,2) <= 0.1);
%
%         marginally_significant_cluster_lh = marginally_significant_cluster_lh(~ismember(marginally_significant_cluster_lh,significant_cluster_lh));
%         marginally_significant_cluster_rh = marginally_significant_cluster_rh(~ismember(marginally_significant_cluster_rh,significant_cluster_rh));
%
%         outputFileName_sig = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName '_significant_cluster']);
%         outputFileName_mar_sig = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName '_mar_significant_cluster']);
%
%         temp_lh = lh_cluster;
%         temp_rh = rh_cluster;
%
%         temp_lh.data(~ismember(temp_lh.data,marginally_significant_cluster_lh)) = 0;
%         temp_lh.data(ismember(temp_lh.data,marginally_significant_cluster_lh)) = 1;
%         temp_rh.data(~ismember(temp_rh.data,marginally_significant_cluster_rh)) = 0;
%         temp_rh.data(ismember(temp_rh.data,marginally_significant_cluster_rh)) = 1;
%
%         temp_vol_lh = lh_Vol;
%         temp_vol_rh = rh_Vol;
%
%         temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
%         temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;
%
%         temp_vol_lh.tstep = tStep ./ 1000;
%         temp_vol_rh.tstep = tStep ./ 1000;
%
%         mne_write_stc_file1([outputFileName_mar_sig, '-lh.stc'], temp_vol_lh);
%         mne_write_stc_file1([outputFileName_mar_sig, '-rh.stc'], temp_vol_rh);
%
%         temp_lh = lh_cluster;
%         temp_rh = rh_cluster;
%
%         temp_lh.data(~ismember(temp_lh.data,significant_cluster_lh)) = 0;
%         temp_lh.data(ismember(temp_lh.data,significant_cluster_lh)) = 1;
%         temp_rh.data(~ismember(temp_rh.data,significant_cluster_rh)) = 0;
%         temp_rh.data(ismember(temp_rh.data,significant_cluster_rh)) = 1;
%
%         temp_vol_lh = lh_Vol;
%         temp_vol_rh = rh_Vol;
%
%         temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
%         temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;
%
%         temp_vol_lh.tstep = tStep ./ 1000;
%         temp_vol_rh.tstep = tStep ./ 1000;
%
%         mne_write_stc_file1([outputFileName_sig, '-lh.stc'], temp_vol_lh);
%         mne_write_stc_file1([outputFileName_sig, '-rh.stc'], temp_vol_rh);
%
%         %% Generate report of the cluster level stats
%         report =cell(1,8);
%         for i = 1:size(significant_cluster_lh,1)
%             report(i,1) = {significant_cluster_lh(i)}; % cluster number
%             report(i,2) = {(find(sum(lh_cluster.data ==significant_cluster_lh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
%             report(i,3) = {(find(sum(lh_cluster.data ==significant_cluster_lh(i))>0, 1, 'last' )-1)*tStep+thisTimeWindow(1)}; % end time in ms
%             report(i,4) = {(find(sum(lh_cluster.data ==significant_cluster_lh(i))==max(sum(lh_cluster.data ==significant_cluster_lh(i))),1,'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
%             temp_lh = lh_cluster;
%             temp_lh.data(~ismember(temp_lh.data,significant_cluster_lh(i))) = 0;
%             temp_lh.data(ismember(temp_lh.data,significant_cluster_lh(i))) = 1;
%             temp_vol_lh = lh_Vol;
%             temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
%             [x,y]=find(temp_vol_lh.data==max(max(temp_vol_lh.data)));
%             report(i,5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
%             report(i,6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
%             report(i,7) = {lh_cluster_size(significant_cluster_lh(i),1)}; % cluster stats, i.e. cluster mass
%             report(i,8) = {lh_cluster_size(significant_cluster_lh(i),2)}; % cluster level p value after the correction
%         end
%         for i = 1:size(marginally_significant_cluster_lh,1)
%             report(i+size(significant_cluster_lh,1),1) = {marginally_significant_cluster_lh(i)}; % cluster number
%             report(i+size(significant_cluster_lh,1),2) = {(find(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
%             report(i+size(significant_cluster_lh,1),3) = {(find(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))>0, 1, 'last' )-1)*tStep+thisTimeWindow(1)}; % end time in ms
%             report(i+size(significant_cluster_lh,1),4) = {(find(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))==max(sum(lh_cluster.data ==marginally_significant_cluster_lh(i))),1, 'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
%             temp_lh = lh_cluster;
%             temp_lh.data(~ismember(temp_lh.data,marginally_significant_cluster_lh(i))) = 0;
%             temp_lh.data(ismember(temp_lh.data,marginally_significant_cluster_lh(i))) = 1;
%             temp_vol_lh = lh_Vol;
%             temp_vol_lh.data = temp_vol_lh.data .* temp_lh.data;
%             [x,y]=find(temp_vol_lh.data==max(max(temp_vol_lh.data)));
%             report(i+size(significant_cluster_lh,1),5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
%             report(i+size(significant_cluster_lh,1),6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
%             report(i+size(significant_cluster_lh,1),7) = {lh_cluster_size(marginally_significant_cluster_lh(i),1)}; % cluster stats, i.e. cluster mass
%             report(i+size(significant_cluster_lh,1),8) = {lh_cluster_size(marginally_significant_cluster_lh(i),2)}; % cluster level p value after the correction
%         end
%
%         % alternate file writing in table format
%         fileID = fopen(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_left-report.txt']),'w');
%         if fileID~=-1 && ~isempty(report)
%             fprintf(fileID, '%5i %30s %36s %36s %44s %55s %50s %50s %55s\n','', 'c1: cluster number','c2: start time in ms', 'c3: end time in ms', 'c4: peak time with maximum space', 'c5: peak vertex number with maximum fit, i.e. r value',...
%                 'c6:peak time with maximum fit, i.e. r value', 'c7: cluster stats, i.e. cluster mass', 'c8: cluster level p value after the correction');
%             totalClusters = size(report);
%             for i = 1:totalClusters(1)
%                 fprintf(fileID, '%5i %25g %36g %36g %44g %55g %50g %50g %55g\n', i, flipud(cell2mat(report(i,:))));
%             end
%         end
%
%         % if ~isempty(report)
%         %     dlmwrite(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_left-report.txt']),cell2mat(report),'\t');
%         % end
%
%         report =cell(1,8);
%         for i = 1:size(significant_cluster_rh,1)
%             report(i,1) = {significant_cluster_rh(i)}; % cluster number
%             report(i,2) = {(find(sum(rh_cluster.data ==significant_cluster_rh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
%             report(i,3) = {(find(sum(rh_cluster.data ==significant_cluster_rh(i))>0, 1, 'last')-1)*tStep+thisTimeWindow(1)}; % end time in ms
%             report(i,4) = {(find(sum(rh_cluster.data ==significant_cluster_rh(i))==max(sum(rh_cluster.data ==significant_cluster_rh(i))),1, 'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
%             temp_rh = rh_cluster;
%             temp_rh.data(~ismember(temp_rh.data,significant_cluster_rh(i))) = 0;
%             temp_rh.data(ismember(temp_rh.data,significant_cluster_rh(i))) = 1;
%             temp_vol_rh = rh_Vol;
%             temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;
%             [x,y]=find(temp_vol_rh.data==max(max(temp_vol_rh.data)));
%             report(i,5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
%             report(i,6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
%             report(i,7) = {rh_cluster_size(significant_cluster_rh(i),1)}; % cluster stats, i.e. cluster mass
%             report(i,8) = {rh_cluster_size(significant_cluster_rh(i),2)}; % cluster level p value after the correction
%         end
%         for i = 1:size(marginally_significant_cluster_rh,1)
%             report(i+size(significant_cluster_rh,1),1) = {marginally_significant_cluster_rh(i)}; % cluster number
%             report(i+size(significant_cluster_rh,1),2) = {(find(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))>0, 1, 'first')-1)*tStep+thisTimeWindow(1)}; % start time in ms
%             report(i+size(significant_cluster_rh,1),3) = {(find(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))>0, 1, 'last')-1)*tStep+thisTimeWindow(1)}; % end time in ms
%             report(i+size(significant_cluster_rh,1),4) = {(find(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))==max(sum(rh_cluster.data ==marginally_significant_cluster_rh(i))),1, 'first')-1)*tStep+thisTimeWindow(1)}; % peak time with maximum space
%             temp_rh = rh_cluster;
%             temp_rh.data(~ismember(temp_rh.data,marginally_significant_cluster_rh(i))) = 0;
%             temp_rh.data(ismember(temp_rh.data,marginally_significant_cluster_rh(i))) = 1;
%             temp_vol_rh = rh_Vol;
%             temp_vol_rh.data = temp_vol_rh.data .* temp_rh.data;
%             [x,y]=find(temp_vol_rh.data==max(max(temp_vol_rh.data)));
%             report(i+size(significant_cluster_rh,1),5) = {x(1,1)}; % peak vertex number with maximum fit, i.e. r value
%             report(i+size(significant_cluster_rh,1),6) = {(y(1,1)-1)*tStep+thisTimeWindow(1)}; % peak time with maximum fit, i.e. r value
%             report(i+size(significant_cluster_rh,1),7) = {rh_cluster_size(marginally_significant_cluster_rh(i),1)}; % cluster stats, i.e. cluster mass
%             report(i+size(significant_cluster_rh,1),8) = {rh_cluster_size(marginally_significant_cluster_rh(i),2)}; % cluster level p value after the correction
%         end
%
%         % alternate file writing in table format
%         fileID = fopen(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_right-report.txt']),'w');
%         if fileID~=-1 && ~isempty(report)
%             fprintf(fileID, '%5i %30s %36s %36s %44s %55s %50s %50s %55s\n','', 'c1: cluster number','c2: start time in ms', 'c3: end time in ms', 'c4: peak time with maximum space', 'c5: peak vertex number with maximum fit, i.e. r value',...
%                 'c6:peak time with maximum fit, i.e. r value', 'c7: cluster stats, i.e. cluster mass', 'c8: cluster level p value after the correction');
%             totalClusters = size(report);
%             for i = 1:totalClusters(1)
%                 fprintf(fileID, '%5i %25g %36g %36g %44g %55g %50g %50g %55g\n', i, flipud(cell2mat(report(i,:))));
%             end
%         end
%
%         % if ~isempty(report)
%         %     dlmwrite(fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_right-report.txt']),cell2mat(report),'\t');
%         % end
%
%         disp(['vertex level threshold = ', num2str(vertex_level_threshold)]);
%         nullDistribution = sort(nullDistribution);
%         cluster_level_threshold = nullDistribution(ceil(size(nullDistribution,2)*(1-0.05)));
%         disp(['cluster level threshold = ', num2str(cluster_level_threshold)]);
%         figure;
%         hist(nullDistribution,0:max(nullDistribution)/100:max(nullDistribution));
%         title('null distribution');
%     end % for: mask
% end % if maskingFlag
