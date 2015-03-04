function nearbySensors = getNearbySensors(sensorSites, adjacencyMatrix, radius)

% Regard the danger/
% Of a recursive function/
% Use at your own risk/
%
% Cai Wingfield 5-2010

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if radius ~= floor(radius)
	error('Radius must be an integer.');
end%if

if radius == 0
	nearbySensors = sensorSites;
elseif radius > 0
	
	nSensors = numel(sensorSites);

	nearbySensors = [];
	for sensorIndex = 1:nSensors
		thisSensor = sensorSites(sensorIndex);
		nearbySensors = [nearbySensors thisSensor adjacencyMatrix(thisSensor, :)];
	end%for:sensorIndex
	
	nearbySensors = nearbySensors(~isnan(nearbySensors));
	nearbySensors = unique(nearbySensors);

	% Finally, call the whole thing recursively! O_o
	nearbySensors = getNearbySensors(nearbySensors, adjacencyMatrix, radius - 1);

else
	error('Radius must be positive.');
end%if
