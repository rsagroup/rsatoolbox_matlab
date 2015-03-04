function closeProgressMonitor(h);

global timesAtFirstCall progressMonitorHandles;
ind=find(progressMonitorHandles==h);
progressMonitorHandles(ind)=[];
timesAtFirstCall(ind,:)=[];
close(h);
