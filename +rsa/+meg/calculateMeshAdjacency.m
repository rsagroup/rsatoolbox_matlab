function searchlightAdjacency = calculateMeshAdjacency(nVertices, searchlightRadius_mm, userOptions)

% searchlightAdjacency = calculateMeshAdjacency(nVertices, searchlightRadius_mm, userOptions)
%
% All credit to Su Li and Andy Thwaites for working out how to do this and writing the original implementation
% CW 5-2010, last updated by Li Su - 1 Feb 2012

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd;

maximumVertices = 40968;
downsampleRate = log(ceil(maximumVertices / nVertices)) / log(4) * 4; 
freesurferResolution = 1.25; % mm between freesurfer vertices

searchlightRadius_freesurfer = searchlightRadius_mm / freesurferResolution;

searchlightCircleRadii_MNE = ceil(downsampleRate*(1:(searchlightRadius_freesurfer/downsampleRate)));

matrixFilename = [userOptions.analysisName '_vertexAdjacencyTable_radius-' num2str(searchlightRadius_mm) 'mm_' num2str(nVertices) '-vertices.mat'];

gotoDir(userOptions.rootPath, 'ImageData');

if ~exist(matrixFilename, 'file')

	fprintf(['The file "' matrixFilename ' doesn''t exist yet, creating it...\n']);
	
	% building a hash table to store adjacency information of all vertexs. 
	hashTableL = findAdjacentVerts(userOptions.averageSurfaceFile); % the resulting hash table is ht_*
	% I might be wrong, but I have discovered that the adjacent
	% vertex indexing is the same across hemispheres, so I only do
	% it for the left. Li Su
	
	fprintf('Building vertex adjacency matrix...\n');
		
	for currentSearchlightCentre = 1:nVertices
	
		if mod(currentSearchlightCentre,floor(nVertices/11)) == 0
			fprintf(['   Working on the vertex ' num2str(currentSearchlightCentre) ' of ' num2str(nVertices) ': ' num2str(floor(100*(currentSearchlightCentre/nVertices))) '%%\n']);
		end
		
		verticesWithinSearchlight = [];
			
		for rMNE = 1:numel(searchlightCircleRadii_MNE)
			freesurferVerticesWithinThisMNERadius = getadjacent(num2str(currentSearchlightCentre),searchlightCircleRadii_MNE(rMNE),hashTableL);
			verticesWithinSearchlight = [verticesWithinSearchlight; freesurferVerticesWithinThisMNERadius(freesurferVerticesWithinThisMNERadius <= nVertices)]; % By removing any which are greater than nVertices, we effectively downsample by the necessary ammount.  This seems a little too clever to work? < or <=?
		end
		
		searchlightAdjacency(currentSearchlightCentre,1:numel(verticesWithinSearchlight)) = verticesWithinSearchlight';
	end
	
	fprintf('      Done!\n');
	
	searchlightAdjacency(searchlightAdjacency == 0) = NaN;
	
	% Save this matrix
	
	cd(fullfile(userOptions.rootPath, 'ImageData'));
	
	save(matrixFilename, 'searchlightAdjacency');

else

	fprintf(['The file "' matrixFilename ' has already been created, loading it...\n']);

	load(matrixFilename);

end

cd(returnHere);
