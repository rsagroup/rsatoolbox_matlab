function timeString=sec2daysHrsMinSec_string(timePeriod_s)

% FUNCTION
%       returns as a string in terms of days, hours, minutes and
%       seconds the time period timePeriod_s expressed in seconds.
%
% USAGE
%       [days,hours,minutes,seconds,secondfraction]=sec2daysHrsMinSec(timePeriod_s)

[days,hours,minutes,seconds,secondfraction]=sec2daysHrsMinSec(timePeriod_s)

if days>1
    days_str=[num2str(days),' days, '];
elseif days==1
    days_str=[num2str(days),' day, '];
else
    days_str=[];
end

if hours>1
    hours_str=[num2str(hours),' hours, '];
elseif hours==1
    hours_str=[num2str(hours),' hour, '];
else
    hours_str=[];
end

if minutes>1
    minutes_str=[num2str(minutes),' minutes, '];
elseif minutes==1
    minutes_str=[num2str(minutes),' minute, '];
else
    minutes_str=[];
end

if seconds>1
    seconds_str=[num2str(seconds),' seconds'];
elseif seconds==1
    seconds_str=[num2str(seconds),' second'];
else
    seconds_str='less than a second';
end

timeString=[days_str,hours_str,minutes_str,seconds_str];

end%function

