function MDSConditions(RDMs, userOptions, localOptions)

% MDSConditions(RDMs, userOptions[, localOptions])
% This function performs first-order MDS, visualising the dissimilarities
% between brain response patterns. If specified in userOptions, a coloured
% dot will be drawn for each condition and a text label is displayed for
% each condition as well. The user can omit the text labels by defining
% them as empty strings. The function will generate a separate MDS plot for
% each RDM in the argument RDMs.
%
% RDMs: A structure of RDMs. All RDMs in here will be concatenated and have
% their condition MDS plots displayed.
% userOptions: The options structure.
%   userOptions.analysisName: A string which is prepended to the saved files.
%   userOptions.rootPath: A string describing the root path where files will
% be saved (inside created directories). userOptions.saveFigurePDF: A
% Boolean value. If true, the figure is saved as a PDF. Defaults to false.
%   userOptions.saveFigurePS: A Boolean value. If true, the figure is saved
% as a PS. Defaults to false. userOptions.saveFigureFig: A Boolean value.
% If true, the figure is saved as a MATLAB .fig file. Defaults to false.
%   userOptions.displayFigures: A Boolean value. If true, the figure remains
% open after it is created. Defaults to true. 
%   userOptions.conditionLabels: A cell array containing the names of the
% conditions in this experiment.
%   userOptions.criterion: The criterion which will be minimised to 
% optimise the MDS arrangement. Defaults to metric stress.
%   userOptions.rubberbands: Boolean value. If true, rubberbands indicating
% MDS distortion are drawn on the MDS plot. Defaults to true. 
%   userOptions.conditionColours: a [nConditions 3]-sized matrix indicating
% an [R G B] triple colour for each condition. Defaults to all black.
%   userOptions.convexHulls: a vector of length equal to the number of 
% conditions. Each entry in the vector corresponds to the same-index
% condition, and the number for this entry represents a category for this
% condition. Convex hulls are drawn around all conditions of the same
% category on the MDS plots, coloured by the first colour in
% userOptions.conditionColours for the points in this category. For
% example: [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4] would represent four
% categories, each with four conditions. If unset, convex hulls will not be
% drawn. 
%   localOptions: Further options specific to this function.
%    localOptions.titleString: If set, this will replace the default title
%    for
% the dendrogram. localOptions.alternativeConditionLabels: A cell array
% containing alternative names for each condition for display on figures.
%    localOptions.figureNumber: If specified AND if only one figure will be
% produced, this will set the figure number of the produced figure.
% Otherwise the figure number will be randomly generated (and probably
% large).

%
%
% Cai Wingfield 3-2010, 5-2010, 6-2010
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
if isfield(localOptions, 'figureNumber') && numel(RDMs) ~= 1, error('MDSConditions:MultipleFigures', 'Can''t use a single specified figure number if there will be more than one figure created in this function run.'); end%if
if ~isfield(userOptions, 'analysisName'), error('MDSConditions:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('MDSConditions:NoRootPath', 'rootPath must be set. See help'); end%if
if ~isfield(userOptions, 'conditionLabels'), error('MDSConditions:NoConditionLabels', 'conditionLabels must be set. See help.'); end%if
userOptions = setIfUnset(userOptions, 'saveFigurePDF', false);
userOptions = setIfUnset(userOptions, 'saveFigurePS', false);
userOptions = setIfUnset(userOptions, 'saveFigureFig', false);
userOptions = setIfUnset(userOptions, 'displayFigures', true);
userOptions = setIfUnset(userOptions, 'criterion', 'metricstress');
userOptions = setIfUnset(userOptions, 'rubberbands', true);
userOptions = setIfUnset(userOptions, 'conditionColours', zeros(numel(userOptions.conditionLabels), 3));

% Pull the RDMs into a 1-d structured array
RDMs = interleaveRDMs(RDMs);

% How many RDMs?
nRDMs = numel(RDMs);

figureNumberBase = 1000000*floor(100*rand);

% Set MDSOptions
if isfield(localOptions, 'alternativeConditionLabels')
	MDSOptions.textLabels = localOptions.alternativeConditionLabels;
else
	MDSOptions.textLabels = userOptions.conditionLabels;
end
MDSOptions.MDSCriterion = userOptions.criterion;
MDSOptions.rubberbandGraphPlot = userOptions.rubberbands;
MDSOptions.dotColours = userOptions.conditionColours;
if isfield(userOptions, 'convexHulls'), MDSOptions.convexHulls = userOptions.convexHulls; end%if

for RDMi = 1:nRDMs

	if isfield(localOptions, 'figureNumber') && numel(RDMs) == 1
		thisFigureNumber = localOptions.figureNumber;
	else
		thisFigureNumber = figureNumberBase + RDMi;
	end%if
	
	veryLocalOptions = MDSOptions;
	veryLocalOptions.figI_textLabels = [thisFigureNumber 2 1 2];
	veryLocalOptions.figI_shepardPlots = [thisFigureNumber 2 1 1];

	RDMName = removeSpaces(strrep(RDMs(RDMi).name, '|', ':'));

	veryLocalOptions.fileName = ['FirstOrderMDS_' RDMName]; % Just what's after "userOptions.analysisName_"

	if isfield(localOptions, 'titleString')
		veryLocalOptions.titleString = localOptions.titleString;
	else
		veryLocalOptions.titleString = ['Conditions MDS for RDM "' RDMs(RDMi).name '"'];
	end%if:localOptions.titleString
	
	fprintf(['Drawing MDS arrangement for conditions of for RDM "' RDMName '"...\n        "' veryLocalOptions.titleString '" [figure ' num2str(thisFigureNumber) ']\n']);

	figureMDSArrangement(RDMs(RDMi), userOptions, veryLocalOptions);

end%for

cd(returnHere);
