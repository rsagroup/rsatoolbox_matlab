function closeProgressMonitors;

global timesAtFirstCall progressMonitorHandles;

for progressMonitorI=1:numel(progressMonitorHandles)
    h=progressMonitorHandles(progressMonitorI);
    if ishandle(h)
        close(h);
    end
end
clear progressMonitorHandles timesAtFirstCall;

