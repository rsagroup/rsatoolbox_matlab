function options=getUserOptions(newOptions,defOptions,allowed);
% function Opt=getUserOptions(options,Opt);
% 
% Deals with user-option structure. It allows the user to specify user
% options as a structure, or as a set of fieldnames and values
%
% INPUTS
%   options: either a structure of user options or a cell array
%            of varargins with 'fieldname1',value1,'fieldname2',value2 ...
% OPTIONAL:
%   Opt:     structure of default options
%   allowed: Cell array of allowed fieldnames for error checking
% 
% EXAMPLE:
% function myFunc(input,varargin);
%   % Specify default options
%   Opt.wdir = 'test';
%   Opt.method = 1;
%   % Allow user to overwrite defaults
%   Opt = rsa_getUserOptions(varargin,Opt,{'wdir','method');
%
% % User can call this with
%   myFunc(input,myOpts) Or
%   myFunc(input,'wdir','root','method',2);
%   myFunc(input,myOpts,'fieldname',value,otherOptions);
% Joern Diedrichsen
% j.diedrichsen@ucl.ac.uk
% 2/2015

if (nargin<3)
    allowed=[];
end;
if (nargin<2)
    defOptions = [];
end;
options = defOptions; 

% Deal with option structure
c=1;
while c<=length(newOptions)
    if (isstruct(newOptions{c}))
        names = fieldnames(newOptions{c});
        for f=1:length(names);
            if ~isempty(allowed)
                if ~any(strcmp(names{f},allowed));
                    msg = sprintf('User argument ''%s'' is not allowed',names{f})
                    msg = [msg sprintf('Allowed names are:\n')];
                    msg = [msg sprintf('%s\n',allowed{:})];
                    error(msg);
                end;
            end;
            options.(names{f})=newOptions{c}.(names{f});
        end;
        c=c+1;
        % Alternatively, deal with number of input strings
    elseif (ischar(newOptions{c}))
        if ~isempty(allowed)
            if ~any(strcmp(newOptions{c},allowed));
                msg = sprintf('User argument ''%s'' is not allowed',newOptions{c})
                msg = [msg sprintf('Allowed names are:\n')];
                msg = [msg sprintf('%s\n',allowed{:})];
                error(msg);
            end;
        end;
        if length(newOptions)==c
            error(sprintf('Option ''%s'' must be followed by a value',newOptions{c}));
        end;
        options.(newOptions{c})=newOptions{c+1};
        c=c+2;
    else
        error('Options must be either fieldnames or structures');
    end;
end;
