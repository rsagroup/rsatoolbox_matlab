function dendrogramConditions(RDMs, userOptions, localOptions)

% Will generate a text-labelled dendrogram for each RDM in RDMs.
%
% dendrogramConditions(RDMs, userOptions[, localOptions])
%
%        RDMs --- The RDMs to be operated on.
%                RDMs is some kind of struct containing RDMs. It can be of
%                whatever dimension but must have fields:
%                        RDM
%                        name
%                        color
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                        userOptions.conditionLabels. Defaults to false.
%                userOptions.conditionColours
%                        A [nConditions 3]-sized matrix indicating an [R G B]
%                        tripple colour for each condition.
%                userOptions.saveFigurePDF
%                        A boolean value. If true, the figure is saved as a PDF.
%                        Defaults to false.
%                userOptions.saveFigurePS
%                        A boolean value. If true, the figure is saved as a PS.
%                        Defaults to false.
%                userOptions.saveFigureFig
%                        A boolean value. If true, the figure is saved as a
%                        MATLAB .fig file. Defaults to false.
%                userOptions.displayFigures
%                        A boolean value. If true, the figure remains open after
%                        it is created. Defaults to true.
%
%        localOptions --- Further options.
%                localOptions.titleString
%                        If set, this will replace the default title for the
%                        dendrogram.
%                localOptions.alternativeConditionLabels
%                        A cell array containing alternative names for each
%                        condition for display on figures. Defaults to be the
%                        same as userOptions.conditionLables.
%                localOptions.useAlternativeConditionLabels
%                        Boolean value. If true,
%                        localOptions.alternativeConditionLabels are used
%                        to label the diagram instead of
%                        userOptions.conditionLabels.
%                localOptions.figureNumber
%                        If specified AND if only one figure will be produced,
%                        this will set the figure number of the produced figure.
%                        Otherwise the figure number will be randomly generated
%                        (and probably large).
%
% Figures may be saved according to preferences.
%
% Cai Wingfield 3-2010, 5-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd;

%% Set defaults and check options struct
if nargin == 2, localOptions = struct(); end%if:nargin
if ~isfield(userOptions, 'analysisName'), error('figureDendrogram:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('figureDendrogram:NoRootPath', 'rootPath must be set. See help'); end%if
localOptions = setIfUnset(localOptions, 'useAlternativeConditionLabels', false);
if localOptions.useAlternativeConditionLabels && ~isfield(localOptions, 'alternativeConditionLabels'), error('dendrogramConditions:NoAlternativeConditionLablels', 'Can''t use alternative condition lables because there aren''t any.'); end%if
if isfield(localOptions, 'figureNumber') && numel(RDMs) ~= 1, error('dendrogramConditions:MultipleFigures', 'Can''t use a single specified figure number if there will be more than one figure created in this function run.'); end%if
userOptions = setIfUnset(userOptions, 'saveFigurePDF', false);
userOptions = setIfUnset(userOptions, 'saveFigurePS', false);
userOptions = setIfUnset(userOptions, 'saveFigureFig', false);
userOptions = setIfUnset(userOptions, 'displayFigures', true);

% Pull the RDMs into a 1-d structured array
RDMs = interleaveRDMs(RDMs);

% How many RDMs?
nRDMs = numel(RDMs);

figureNumberBase = 1000000*floor(100*rand);

for RDMi = 1:nRDMs

	thisFigureNumberBase = figureNumberBase + RDMi*1000;

	thisFigureNumber = thisFigureNumberBase + 2;

	if isfield(localOptions, 'figureNumber') && numel(RDMs) == 1
		veryLocalOptions.figureNumber = localOptions.figureNumber;
	else
		veryLocalOptions.figureNumber = thisFigureNumber;
	end%if

	RDMName = removeSpaces(strrep(RDMs(RDMi).name, '|', ':'));

	veryLocalOptions.fileName = ['Dendrogram_' RDMName]; % Just what's after "userOptions.analysisName_"
	veryLocalOptions.linkageType = 'average';
	if localOptions.useAlternativeConditionLabels
		veryLocalOptions.labels = localOptions.alternativeConditionLabels;
	else
		veryLocalOptions.labels = userOptions.conditionLabels;
	end
	
	if isfield(localOptions, 'titleString')
		veryLocalOptions.titleString = localOptions.titleString;
	else
		veryLocalOptions.titleString = ['Conditions dendrogram for RDM "' RDMs(RDMi).name '"'];
	end%if:localOptions.titleString

	fprintf(['Drawing dendrogram for RDM "' RDMName '"...\n        "' veryLocalOptions.titleString '" [figure ' num2str(thisFigureNumber) ']\n']);

	figureDendrogram(RDMs(RDMi).RDM, userOptions, veryLocalOptions);

end%for

cd(returnHere);
