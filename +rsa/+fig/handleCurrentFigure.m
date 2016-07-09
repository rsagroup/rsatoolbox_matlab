function handleCurrentFigure(fileName, userOptions)

% handleCurrentFigure(fileName, userOptions)
% 
% handleCurrentFigure is a function which will, based on the preferences set in
% userOptions, save the current figure as a .pdf, a .ps or a .fig; either
% leaving it open or closing it.
%
%        fileName --- The name of the file to be saved.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.saveFiguresPDF
%                        A boolean value. If true, the figure is saved as a PDF.
%                        Defaults to false.
%                userOptions.saveFiguresPS
%                        A boolean value. If true, the figure is saved as a PS.
%                        Defaults to false.
%                userOptions.saveFiguresFig
%                        A boolean value. If true, the figure is saved as a
%                        MATLAB .fig file. Defaults to false.
%                userOptions.saveFiguresEps
%                        A boolean value. If true, the figure is saved as an
%                        EPS. Defaults to false.
%                userOptions.saveFiguresJpg
%                        A boolean value. If true, the figure is saved as a jpeg
%                        image. (Thanks to Alex for the tip!)
%                userOptions.displayFigures
%                        A boolean value. If true, the figure remains open after
%                        it is created. Defaults to true.
%
% Cai Wingfield 3-2010
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

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), userOptions.analysisName='unnamed analysis'; end%if
if ~isfield(userOptions, 'rootPath'), userOptions.rootPath=pwd; end%if
userOptions = setIfUnset(userOptions, 'saveFiguresPDF', false);
userOptions = setIfUnset(userOptions, 'saveFiguresPS', false);
userOptions = setIfUnset(userOptions, 'saveFiguresFig', false);
userOptions = setIfUnset(userOptions, 'saveFiguresEps', false);
userOptions = setIfUnset(userOptions, 'saveFiguresJpg', false);
userOptions = setIfUnset(userOptions, 'displayFigures', true);

appendFlag = 0;

if userOptions.saveFiguresPDF
	exportCurrentFigAsPDF(fileName, userOptions);
end%if:pdf
if userOptions.saveFiguresPS
	exportCurrentFigAsPostscript(fileName, appendFlag, userOptions);
end%if:ps
if userOptions.saveFiguresFig
	hgsave(fileName);
end%if:fig
if userOptions.saveFiguresEps
	exportfig(gcf, fileName, 'Format', 'eps', 'Color', 'cmyk', 'Renderer', 'zbuffer', 'Resolution', 300);
end%if:eps
if userOptions.saveFiguresJpg
	saveas(gcf, [fileName '.jpg']);
end%if:jpg
if ~userOptions.displayFigures
	close;
end%if:display

end%function
