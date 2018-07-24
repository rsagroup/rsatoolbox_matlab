function figureRDMs(RDMs, userOptions, localOptions)
%  figureRDMs(RDMs, userOptions[, localOptions])
%
%  figureInterleavedRDMs is a function which accepts a multidimensional struct of
%  RDMs, puts them into a 1D struct, and then shows them.
%
%        RDMs --- A struct of RDMs.
%                All RDMs in here will be concatenated and displayed.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
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
%                userOptions.imagelables
%                        Defaults to empty (no image labels).
%                userOptions.rankTransform
%                        Boolean value. If true, values in each RDM are
%                        separately transformed to lie uniformly in [0,1] with
%                        ranks preserved. If false, true values are represented.
%                        Defaults to true. This is for display only.
%                userOptions.colourScheme
%                        A colour scheme for the RDMs. For Matlab version
%                        before R2015b,
%                        defualts to jet(64); for versions afterward,
%                        defaults to colormap(jet(256)).
%                        
%        localOptions --- Further options.
%                localOptions.figureNumber
%                        If specified, this will set the figure number of the
%                        produced figure. Otherwise the figure number will be
%                        randomly generated (and probably large).
%                localOptions.fileName
%
% May save figures according to preferences.
%
%  Cai Wingfield 11-2009, 6-2010
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
% Detect Matlab version (R2015b changes default color maps: should use
% colormap() to specify the color of the color bar.
matlabVersion = str2num(cell2mat(regexp(version('-release'),'\d+','match')));
if matlabVersion>2015
    color=colormap(jet(256));
else
    color=jet(64);
end
%% Set defaults and check options struct
if nargin == 2, localOptions = struct(); end%if:nargin
if ~isfield(userOptions, 'analysisName'), error('figureInterleavedRDMs:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('figureInterleavedRDMs:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'saveFigurePDF', false);
userOptions = setIfUnset(userOptions, 'saveFigurePS', false);
userOptions = setIfUnset(userOptions, 'saveFigureFig', false);
userOptions = setIfUnset(userOptions, 'displayFigures', true);
userOptions = setIfUnset(userOptions, 'rankTransform', true);
userOptions = setIfUnset(userOptions, 'colourScheme', color);

if ~isfield(userOptions,'dpi'), userOptions.dpi=300; end;
if ~isfield(userOptions,'tightInset')
    userOptions.tightInset = false;
end

appendFlag = 0; % Put this into a prompt depending on existant files?

if ~isfield(localOptions, 'figureNumber')
	localOptions.figureNumber = 1;
end
if ~isfield(localOptions, 'fileName')
	localOptions.fileName = 'RDMs';
end

RDMs = interleaveRDMs(RDMs); % Pull the RDMs into a 1-d structured array

%% Now display

% This is not available on mac, it causes an error in Matlab 2014
% opengl software

if isfield(userOptions, 'imagelabels')
    imagelabels = userOptions.imagelabels;
else
    imagelabels = [];
end

if userOptions.rankTransform
	showRDMs(RDMs, localOptions.figureNumber, true, [0 1], true, [], imagelabels, userOptions.colourScheme);
else
	showRDMs(RDMs, localOptions.figureNumber, false, [], true, [], imagelabels, userOptions.colourScheme);
end%rankTransform

gotoDir(userOptions.rootPath, 'Figures');
fileName = [userOptions.analysisName '_' localOptions.fileName];
handleCurrentFigure(fileName, userOptions);

cd(returnHere);
