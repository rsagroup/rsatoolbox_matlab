function betas = getDataFromSPM(userOptions)
%
% getDataFromSPM is a function which will extract from the SPM metadata the
% correspondence between the beta image filenames and the condition and session
% number.
%
% function betas = getDataFromSPM(userOptions)
%
%        betas --- The array of info.
%               betas(condition, session).identifier is a string which referrs
%               to the filename (not including path) of the SPM beta image.
%
%        userOptions --- The options struct.
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names.
%                userOptions.betaPath
%                        A string which contains the absolute path to the
%                        location of the beta images. It can contain the
%                        following wildcards which would be replaced as
%                        indicated:
%                                [[subjectName]]
%                                        To be replaced with the contents of
%                                        subject.
%                                [[betaIdentifier]]
%                                        To be replaced by filenames returned in
%                                        betas.
%                userOptions.conditionLabels
%                        A cell array containing the names of the conditions in
%                        this experiment. Here, these will be used to find
%                        condition response beta predictor images (so make sure
%                        they're the same as the ones used for SPM!).
%  
%  Cai Wingfield 12-2009, 6-2010, 8-2010

	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

	%% Set defaults and check for problems.
	if ~isfield(userOptions, 'betaPath'), error('getDataFromSPM:NoBetaPath', 'userOptions.betaPath is not set. See help.'); end%if
	if ~isfield(userOptions, 'subjectNames'), error('getDataFromSPM:NoSubjectNames', 'userOptions.subjectNames is not set. See help.'); end%if
	if ~isfield(userOptions, 'conditionLabels'), error('getDataFromSPM:NoConditionLabels', 'userOptions.conditionLabels is not set. See help.'); end%if

    firstSubject = userOptions.subjectNames{1};

	readFile = replaceWildcards(fullfile(userOptions.betaPath, 'SPM.mat'), '[[subjectName]]', firstSubject, '[[betaIdentifier]]', '');
	load(readFile);
	nBetas = max(size(SPM.Vbeta));

	nConditions = numel(userOptions.conditionLabels);

	% Extract all info

	highestSessionNumber = 0;

	for betaNumber = 1:nBetas % For each beta file...

		% Get the description of the beta file
		thisBetaDescrip = SPM.Vbeta(betaNumber).descrip;

		% Extract from it the file name, the session number and the condition name
		[thisBetaName, thisSessionNumber, thisConditionName] = extractSingleBetaInfo(thisBetaDescrip);

		% Check if it's one of the conditions
		thisBetaIsACondition = false; % (What a delicious variable name!)
		for condition = 1:nConditions

			thisBetaIsACondition = strcmpi(thisConditionName, userOptions.conditionLabels{condition}); % Check if it's a condition specified by the user
			conditionNumber = condition; % Record the condition number that it is
			if thisBetaIsACondition, break; end%if

		end%for

		% Now only proceed if we're looking at a condition
		if thisBetaIsACondition
		
			% Keep a tally of the higest run yet
			highestSessionNumber = max(highestSessionNumber, thisSessionNumber);
			
			% Store the file name in the betas struct
			betas(thisSessionNumber, conditionNumber).identifier = [thisBetaName '.img'];
			
		end%if
	end%for

end%function

%% Subfunctions: %%

% spm_spm:beta (0001) - Sn(1) all_events*bf(1)
function [betaName, sessionNumber, conditionName] = extractSingleBetaInfo(strIn)

	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

	openBrackets = findstr(strIn, '(');
	closedBrackets = findstr(strIn, ')');
	star = findstr(strIn, '*');

	inFirstBrackets = strIn(openBrackets(1)+1:closedBrackets(1)-1);
	betaName = ['beta_' inFirstBrackets];

	inSecondBrackets = strIn(openBrackets(2)+1:closedBrackets(2)-1);
	sessionNumber = str2num(inSecondBrackets);

	conditionName = strIn(closedBrackets(2)+2:star-1);
end%function
