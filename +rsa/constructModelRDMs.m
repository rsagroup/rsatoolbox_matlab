function [varargout] = constructModelRDMs(rawModels, userOptions)
%
%  constructModelRDMs is a function which parses the User/modelRDMs.m file and
%  saves a struct which is readable by the body other toolbox functions.
%
% [RDMs =] constructModelRDMs(rawModels, userOptions)
%
%        rawModels --- Naturally specified models.
%                A struct in which rawModels.(modelName) is the model RDM.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.ModelColor
%                        A triple indicating the [R G B] value of the colour
%                        which should be used to indicated model RDMs on various
%                        diagrams. Defaults to black ([0 0 0]).
%
% The following files are saved by this function.
%        userOptions.rootPath/RDMs/
%                userOptions.analysisName_Models.mat
%                        Contains a structure of model RDMs, one for each of the
%                        ones from rawModels, but in the form preferred by the
%                        toolbox.
%  
%  Cai Wingfield 11-2009
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd;

if ~isfield(userOptions, 'analysisName'), error('constructModelRDMs:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('constructModelRDMs:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'ModelColor', [0 0 0]);

ModelsFilename = [userOptions.analysisName, '_Models.mat'];

fprintf('Constructing models based on modelRDMs.m...\n');

% Parse
modelNames = fields(rawModels);
nModelNames = numel(modelNames);

% Reconstruct
for model = 1:nModelNames
	Models(model).name = underscoresToSpaces(modelNames{model});
	Models(model).RDM = rawModels.(modelNames{model});
	Models(model).color = userOptions.ModelColor;
end%for

% And save
% fprintf(['Saving Model RDMs to ' fullfile(userOptions.rootPath, 'RDMs', ModelsFilename) '\n']);
disp(['Saving Model RDMs to ' fullfile(userOptions.rootPath, 'RDMs', ModelsFilename)]);
gotoDir(userOptions.rootPath, 'RDMs');
save(ModelsFilename, 'Models');

if nargout == 1
	varargout{1} = Models;
elseif nargout > 0
	error('0 or 1 arguments out, please.');
end%if:nargout

cd(returnHere);

end%function
