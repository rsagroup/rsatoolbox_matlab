% function [stamped_message =] prints(FORMAT, A, ...)
%
% Sends the specified message to the output window, preceeded by a
% timestamp, and followed by a newline.  Optionally returns the string,
% minus the newline.
%
% Treats its arguments just as sprintf does.
%
% Cai Wingfield 2015-03

function stamped_message = prints(varargin)

    % Get the current time asap.
    datestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
    
    % Apply the formatting as supplied in the argumets.
    message = sprintf(varargin{:});
    
    % Build the stamped string.
    stamped_message = ['[', datestamp, '] ', message];
    
    % Send it to the output window, including a trailing newline.
    fprintf([stamped_message, '\n']);
    
end%function
