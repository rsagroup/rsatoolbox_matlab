function D=struct2dataframe(S);
% function D=struct2dataframe(S);
% Converts a structure-array into a dataframe structure.
% The pre-requisite is that all fields in the structure a 1xK matrices or
% cell arrays, so that they can be vertically concatenated
fields = fieldnames(S(1));
numStruct = length(S); 
for f=1:length(fields)
    if (ischar(S(1).(fields{f})))
        for s=1:numStruct
            D.(fields{f}){s,1} = S(s).(fields{f});
        end;
    else
        D.(fields{f}) = vertcat(S(:).(fields{f}));
    end;
    if size(D.(fields{f}))~=numStruct
        warning(sprintf('Ignoring field %s: needs to have 1 row per entry',fields{f})); 
        D = rmfield(D,fields{f}); 
    end; 
end;
