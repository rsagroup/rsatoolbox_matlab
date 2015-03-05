function [days,hours,minutes,seconds,secondfraction]=sec2daysHrsMinSec(timePeriod_s)

% FUNCTION
%       returns and optionally outputs in terms of days, hours, minutes and
%       seconds the time period timePeriod_s expressed in seconds.
%
% USAGE
%       [days,hours,minutes,seconds,secondfraction]=sec2daysHrsMinSec(timePeriod_s)

seconds=floor(timePeriod_s);
secondfraction=timePeriod_s-seconds;

minutes=floor(seconds/60);
seconds=seconds-minutes*60;

hours=floor(minutes/60);
minutes=minutes-hours*60;

days=floor(hours/24);
hours=hours-days*24;

end%function

