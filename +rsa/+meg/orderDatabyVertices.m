
% function sorts the data elements by vertex list to an ascending order
% vertex list. Added this as lower resolution data wasnt ordered as
% presumed. IZ 06/13

function [ordered_data] = orderDatabyVertices(data, vertices)

if ~issorted(vertices)
    for i=1:max(vertices)
        ordered_data(i,:) = data(vertices==i,:);
    end
else
    ordered_data = data;
end

end