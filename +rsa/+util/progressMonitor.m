function h=progressMonitor(counters, maxCounts, title, h)

% FUNCTION
%       displays waitbar and estimates the time needed to finish a task.
%
% ARGUMENTS
% counters
%       a vector of counters of nested loops ordered from outermost to
%       innermost loop. counters are assumed to start at 1.
%
% maxCounts
%       the maximum values of the counters in corresponding order
%
% title
%       the waitbar title
%
% USAGE
%       h=progressMonitor(ones(size(maxCounts)),maxCounts,title)
%       ...
%       for counter1=1:...
%           for counter2=1:...
%                ...
%                ...
%`              for counterN=1:...
%                   progressMonitor([counter1 counter2...counterN],[maxCount1 maxCount2...maxCountN],title,h)
%                   ...
%                   ...
%               end
%           end
%       end
%       closeProgressMonitor(h);



%% preparations
global timesAtFirstCall progressMonitorHandles;

% columnize counters and maxCounts
if size(counters,2)>size(counters,1)
    counters=counters';
end

if size(maxCounts,2)>size(maxCounts,1)
    maxCounts=maxCounts';
end

n=length(counters);

%% compute the progress [0,1]

counters=counters-1; % conservative: current step not assumed to be done
progress=0;

for i=n:-1:1
    progress=(progress+counters(i))/maxCounts(i);
end

if exist('h','var')  % update existing progress monitor
    
    if progress>0
        cI=find(progressMonitorHandles==h); % this must be unique
        
        elapsedTime_s=etime(clock,timesAtFirstCall(cI,:));
        totalTimeEstimate_s=elapsedTime_s/progress;
        remainingTimeEstimate_s=totalTimeEstimate_s-elapsedTime_s;
 
        if cI==1
            % only display output in command window for the outermost progress monitor...
            disp('___________________________________________________');
            disp(title);
            disp('estimated remaining time needed:');
            sec2daysHrsMinSec(remainingTimeEstimate_s,1);
            disp(datestr(clock));
        else
            fprintf('.');
        end

        % separate waitbars for all nested progress monitors...
        remainingTimeString=sec2daysHrsMinSec_string(remainingTimeEstimate_s);
        rtVec = sec2daysHrsMinSec(remainingTimeEstimate_s);
        remainingTimeString=sprintf()
        newTitle=[title,' (',remainingTimeString,' left.) '];
        % newTitle=[title,' (',remainingTimeString,' left.) ',datestr(clock)];
        %         [progress,h,newTitle]
        waitbar(progress,h,newTitle);
    end

    
else   % new progress monitor
    h=waitbar(0,title);
    set(h,'Resize','on');
    %     rect=get(h,'Position'); rect(3)=rect(3)*2; % widen the waitbar box
    %     set(h,'Position',rect);
    
    if sum(counters~=zeros(size(counters)))
        error('ERROR: All counters should be one upon the initial call.');
    end
    
    if any(find(progressMonitorHandles==h))
        % handle already existed previously, so use that slot to avoid
        % ambiguity
        timesAtFirstCall(progressMonitorHandles==h,:)=clock;
    else
        
        progressMonitorHandles=[progressMonitorHandles;h];
        timesAtFirstCall=[timesAtFirstCall;clock];
    end
    % call closeProgressMonitor(h) to avoid memory leakage
    % (which should be negligible for all practical purposes).
end

