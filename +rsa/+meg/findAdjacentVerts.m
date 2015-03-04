% this function returns a hash table containing the adjacent vertexes for
% each vertex in the brain based on freesurfer cortical mash.
% 
% created by Li Su and Andy Thwaites, last updated by Li Su 01 Feb 2012

function ht = findAdjacentVerts(path)

% addpath /opt/mne/matlab/toolbox/ % CW: path doesn't exist.

[verts,faces] = mne_read_surface(path);
numberOfVerts = max(faces);
numberOfFaces = size(faces);


ht = java.util.Hashtable;

facesduplicate = zeros(length(faces)*3, 3);

for i = 1:length(faces)
    q = length(faces);
    % disp(num2str(i));
    facesduplicate(i,1:3) = [faces(i,1) faces(i,2) faces(i,3)];
    facesduplicate(i+q,1:3) = [faces(i,2) faces(i,1) faces(i,3)];
    facesduplicate(i+(q*2),1:3) = [faces(i,3) faces(i,2) faces(i,1)];
end

sortedfaces = sortrows(facesduplicate,1);

thisface = 1;
adjacent = [];
for i = 1:length(sortedfaces)
    % disp(num2str(i));
    face = sortedfaces(i,1);
    if  (face == thisface)
        key = num2str(face);
        adjacent = [adjacent sortedfaces(i,2)];
        adjacent = [adjacent sortedfaces(i,3)];
    else
        unad = unique(adjacent);
        ht.put(key,unad);
        adjacent = [];
        thisface = face;

        
        % now continue as normall
        key = num2str(face);
        adjacent = [adjacent sortedfaces(i,2)];
        adjacent = [adjacent sortedfaces(i,3)];
    end
end


