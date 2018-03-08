function S=dataframe2struct(D);
% function D=struct2dataframe(S);
% Converts a structure-array into a dataframe structure.
% The pre-requisite is that all fields in the structure a 1xK matrices or
% cell arrays, so that they can be vertically concatenated
fields = fieldnames(D);
numStruct = length(D.(fields{s})); 
for f=1:length(fields)
    if (iscell(D.(fields{f})))
        [S(1:numStruct).(fields{f})]=deal(D.(fields{s}){:}); 
    else 
        for s=1:numStruct
            S(s).(fields{f})=D.(fields{f})(s,:); 
        end;
    end; 
end;
