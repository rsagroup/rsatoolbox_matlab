function MEG_timeseries(Models,indexMasks,userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*


modelNumber = userOptions.modelNumber; 
modelName = spacesToUnderscores(Models(modelNumber).name);
inputFileName = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName '_' modelName '_significant_cluster']);

lh = mne_read_stc_file([inputFileName,'-lh.stc']);
rh = mne_read_stc_file([inputFileName,'-rh.stc']);

vertex_level_threshold = min(min(lh.data(lh.data>0)),min(rh.data(rh.data>0)));
[numberOfVertex, numberOfTimeWindows] = size(lh.data);

nMasks = size(userOptions.maskNames,2);

time_series = -200:5:-200+(numberOfTimeWindows-1)*5;
time_series = [time_series;mean(lh.data);mean(rh.data)];
    
for m = 1:2:nMasks    
    masked_lh.data = MEG_masking(lh.data,indexMasks,m);
    masked_rh.data = MEG_masking(rh.data,indexMasks,m+1);

    time_series_lh = mean(masked_lh.data);
    time_series_rh = mean(masked_rh.data);
    
%     for t = 1:numberOfTimeWindows
%         l = masked_lh.data(:,t);
%         r = masked_rh.data(:,t);
%         time_series_lh(t) = mean(l);
%         time_series_rh(t) = mean(r);

    %     if numel(l(l~=0)) > 2
    %         time_series_lh(t) = mean(l(l~=0));
    %     else time_series_lh(t) = NaN;
    %     end
    % 
    %     if numel(r(r~=0)) > 2
    %         time_series_rh(t) = mean(r(r~=0));
    %     else time_series_rh(t) = NaN;
    %     end
%     end

    time_series_lh(isnan(time_series_lh)) = vertex_level_threshold;
    time_series_rh(isnan(time_series_rh)) = vertex_level_threshold;

    time_series = [time_series;time_series_lh;time_series_rh];
end

    masked_lh.data = MEG_masking(lh.data,indexMasks,21) + MEG_masking(lh.data,indexMasks,23) +MEG_masking(lh.data,indexMasks,25) ;
    masked_rh.data = MEG_masking(rh.data,indexMasks,22) + MEG_masking(rh.data,indexMasks,24) +MEG_masking(rh.data,indexMasks,26) ;

    time_series_lh = mean(masked_lh.data);
    time_series_rh = mean(masked_rh.data);
    
%     for t = 1:numberOfTimeWindows
%         l = masked_lh.data(:,t);
%         r = masked_rh.data(:,t);
%         time_series_lh(t) = mean(l);
%         time_series_rh(t) = mean(r);

    %     if numel(l(l~=0)) > 2
    %         time_series_lh(t) = mean(l(l~=0));
    %     else time_series_lh(t) = NaN;
    %     end
    % 
    %     if numel(r(r~=0)) > 2
    %         time_series_rh(t) = mean(r(r~=0));
    %     else time_series_rh(t) = NaN;
    %     end
%     end

    time_series_lh(isnan(time_series_lh)) = vertex_level_threshold;
    time_series_rh(isnan(time_series_rh)) = vertex_level_threshold;

    time_series = [time_series;time_series_lh;time_series_rh];

outputFileName = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName, '_sig_cluster_time_series.csv']);

csvwrite(outputFileName,time_series);
    
