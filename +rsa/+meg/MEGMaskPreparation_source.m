% MEGMaskPreparation_source will load vertex masks from a directory specified in 
% userOptions and will save a struct of these masks coupled with the time-windows 
% of interest also specified in userOptions.
%
% [indexMasks =] MEGMaskPreparation_source(userOptions)
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.maskPath
%                        A string describing the path to the location of the
%                        files for the definition of RoI masks.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Names can be repeated so that each one is
%                        paired with an entry of userOptions.maskTimeWindows.
%                userOptions.maskTimeWindows
%                        A cell array, the same length as
%                        userOptions.maskNames containing vectors of length 2
%                        containing the start and end of the time window of
%                        interest for the corresponding mask from
%                        userOptions.maskNames.
%
% Cai Wingfield 9-2010, updated by Li Su 3-2012

function [varargout] = MEGMaskPreparation_source(userOptions)

import rsa.*
import rsa.fig.*
import rsa.meg.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

% The analysisName will be used to label the files which are eventually saved.
MasksFilename = [userOptions.analysisName, '_Masks.mat'];
DetailsFilename = [userOptions.analysisName, '_MEGMaskPreparation_source_Details.mat'];

promptOptions.functionCaller = 'MEGMaskPreparation_source';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'ImageData', MasksFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

if userOptions.slidingTimeWindow, overwriteFlag = 1;
else overwriteFlag = overwritePrompt(userOptions, promptOptions); end

if overwriteFlag

	% Data
	nMasks = numel(userOptions.maskNames);
	nTimeWindows = numel(userOptions.maskTimeWindows);
	
	% Test to see if things are broken already
	if nMasks ~= nTimeWindows 
		warning('Not same number of masks and time windows. RoI(fixed and sliding time window) require same number of masks and time windows');
	end%if
	
	%% Get Data
	
	fprintf('Building spatiotemporal masks');

    for maskNumber = 1:nMasks % For each mask...

		% Figure out which mask we're looking at
		mask = userOptions.maskNames{maskNumber};
        if nMasks == nTimeWindows 
            timeWindow = userOptions.maskTimeWindows{maskNumber};
        else
            timeWindow = userOptions.dataPointsSearchlightLimits;
        end
		
		% Load the mask
		readPath = [replaceWildcards(userOptions.maskPath, '[[maskName]]', mask) '.label'];
        label = mne_read_label_file(readPath);
% 		vertexMaskCell = readFileToCell(readPath);
% 		
% 		% Remove the first two (junk) lines
% 		vertexMaskCell = vertexMaskCell(3:end);
% 		
% 		% Get the first column only
% 		for i = 1:numel(vertexMaskCell)
% 			vertexMaskCell_split(i,:) = splitStringAtCharacter(vertexMaskCell{i}, ' ');
% 		end%for:i
% 		vertexMaskCell = vertexMaskCell_split(:, 1);
%         for i = 1:numel(vertexMaskCell)
% 			vertexMaskCell{i} = str2num(vertexMaskCell{i});
% 		end%for:i
% 		
% 		% Convert to column vector
% 		vertexMask = cell2mat(vertexMaskCell);
% 		vertexMask = vertexMask(:);
		
		% Determine if it's left or right
		suffix = mask(end-1:end);
        mask_noDash = dashToUnderscores(mask); % updated by Li Su 2-2012
        
		if strcmpi(suffix, 'lh')
			chi = 'L';
		elseif strcmpi(suffix, 'rh')
			chi = 'R';
		else
			error('MEGMaskPreparation:notLhOrRh', ['The mask ' mask ' does not end in "lh" or "rh", and as such can''t allocated to either the left- or right-hand side of the source mesh!']);
		end%if:suffix
		
		% Store in a struct
		%thisFieldname = dashToUnderscores(['tw_' num2str(timeWindow(1)) '_' num2str(timeWindow(2))]);
		indexMasks.(mask_noDash).maskIndices = label.vertices + 1; % update IZ 02/12 previously: vertexMask + 1;
		indexMasks.(mask_noDash).timeIndices = timeWindow(1):timeWindow(2);
		indexMasks.(mask_noDash).chirality = chi;
		
		fprintf('.');

	end%for:maskNumber
	
	fprintf(':\n');

	%% Save relevant info

	timeStamp = datestr(now);

	fprintf(['Saving masks to ' fullfile(userOptions.rootPath, 'ImageData', MasksFilename) '\n']);
	gotoDir(userOptions.rootPath, 'ImageData');
	save(MasksFilename, 'indexMasks');
	
	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');
	
else
	fprintf(['Loading previously saved masks from ' fullfile(userOptions.rootPath, 'ImageData', MasksFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'ImageData', MasksFilename));
end%if:overwriteFlag

if nargout == 1
	varargout{1} = indexMasks;
elseif nargout > 0
	error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started
