function [RDMs,nRDMs]=unwrapRDMs(RDMs_struct)
% unwraps dissimilarity matrices of a structured array (wrapped RDMs) with
% meta data by extracting the dissimilarity matrices (in square or lower
% triangle form) and lining them up along the third dimension. (if they are
% already in that format they are handed back unchanged.) it must be noted
% that different dissimilarity sets (i.e. lower-triangular or squareform
% RDMs) must be different fields of the input structure. The only
% requirement is that the dissimilarities should be stored in a field
% called 'RDM'.
%
% nRDMs is the number of wrapped RDMs in the input. in case where all the
% input RDMs are in lower-triangular or square form, the output would be in
% the same format as well (RDMs will be stacked along a third dimension).
% if the input RDMs are of mixed type, they're all made square. Also, if
% RDMs in the input are of different sizes, they'll all be returned as
% squares with nans outside. 
%
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if isstruct(RDMs_struct)
    % in struct form
    nRDMs=size(RDMs_struct,2);
	
	mixedTypes = false;
	
	for RDMi = 1:nRDMs
		
		thisRDM = RDMs_struct(RDMi).RDM;
		
		% What type of RDM is this?
		if max(size(thisRDM)) == numel(thisRDM)
			% Then it's in ltv form
			thisType = 'ltv';
			thisSize = (1+sqrt(1-4*numel(thisRDM)))/2;
		else
			thisType = 'sq';
			thisSize = size(thisRDM,1);
		end%if:utv
		
		if RDMi == 1
			% What's the first type?
			firstType = thisType;
			firstSize = thisSize;
			maxSize = thisSize;
		elseif ~strcmp(thisType, firstType) || thisSize ~= firstSize;
			% Is each type the same as the first?
			mixedTypes = true;
			if thisSize ~= firstSize;
				maxSize = max(maxSize, thisSize);
			end%if:not the same size
		end%if:RDMi==1
	
		if ~mixedTypes
			switch thisType
				case 'ltv'
					RDMs(1,:,RDMi) = double(thisRDM);
				case 'sq'
					RDMs(:,:,RDMi) = double(thisRDM);
			end%switch:thisType
		end%if:~mixedTypes
		
	end%for:RDMi
	
	% If previous loop was broken...
	
	if mixedTypes
		clear RDMs;
		% All must be square
		RDMs = nan(maxSize, maxSize, nRDMs);
		for RDMi = 1:nRDMs
			thisRDM = RDMs_struct(RDMi).RDM;
			% What type of RDM is this?
			if max(size(thisRDM)) == numel(thisRDM)
				% Then it's in ltv form
				thisRDM = squareform(thisRDM);
			end%if:utv
			RDMs(1:size(thisRDM,1), 1:size(thisRDM,2), RDMi) = double(thisRDM);
		end%for:RDMi
	end%if:mixedTypes
	
else
    % bare already
    RDMs=RDMs_struct;
    nRDMs=size(RDMs,3);
end

end%function
