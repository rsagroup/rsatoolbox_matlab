function MEG_figures(Models,mask,userOptions)

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

[numberOfVertex, numberOfTimeWindows] = size(lh.data);

if size(mask)>0
    lh.data = MEG_masking(lh.data,mask,1);
    rh.data = MEG_masking(rh.data,mask,-1);
end

timeStep = 50 / userOptions.temporalSearchlightResolution;
temp_lh.tmin = lh.tmin;
temp_rh.tmin = rh.tmin;
temp_lh.tstep = lh.tstep*timeStep;
temp_rh.tstep = rh.tstep*timeStep;
temp_lh.vertices = lh.vertices;
temp_rh.vertices = rh.vertices;

%     for i = 1:numberOfVertex
%         l = lh.data(i,:);
%         r = rh.data(i,:);
%         if numel(l(l~=0)) >1
%             temp_lh.data(i,1) = mean(l(l~=0),2);
%         end
%         if numel(r(r~=0)) >1
%             temp_rh.data(i,1) = mean(r(r~=0),2);
%         end
%     end


x = 1;
for t = 1:timeStep:numberOfTimeWindows
    try
        temp_lh.data(:,x) = mean(lh.data(:,t:t+timeStep),2);
        temp_rh.data(:,x) = mean(rh.data(:,t:t+timeStep),2);
    catch
        temp_lh.data(:,x) = mean(lh.data(:,t:end),2);
        temp_rh.data(:,x) = mean(rh.data(:,t:end),2);
    end

    %temp_lh.data(isnan(temp_lh.data)) = 0;
    %temp_rh.data(isnan(temp_rh.data)) = 0;
    x = x + 1;
    
end

temp_lh.data(:,x+1) = mean(lh.data,2);
temp_rh.data(:,x+1) = mean(rh.data,2);
        
outputFileName = fullfile(userOptions.rootPath, 'Results', [userOptions.analysisName, '_', modelName, '_sig_cluster_50ms_average']);

mne_write_stc_file([outputFileName,'-lh.stc'],temp_lh);
mne_write_stc_file([outputFileName,'-rh.stc'],temp_rh);
