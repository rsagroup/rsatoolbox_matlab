% testSignificance is a function which accepts a struct of RoI RDMs and a struct
% of model RDMs and pairwise tests the significance of their relatedness (as
% calculated by compareRDMs).
%
% testSignificance( ...
%                 {referenceRDMs[, referenceRDMs2, ...]}, ...
%                 {comparedRDMs[, comparedRDMs2, ...]}, ...
%                 userOptions ...
%                )
%
%        referenceRDMs, referenceRDMs2, ... --- The first set of RDMs.
%                Each of these RDMs is some kind of struct containing RDMs. They
%                can be of whatever dimension but must have fields:
%                        RDM
%                        name
%                        color
%
%        comparedRDMs, comparedRDMs2, ... --- The second set of RDMs.
%                Each of these RDMs is some kind of struct containing RDMs. They
%                can be of whatever dimension but must have fields:
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
%                userOptions.significanceTestPermutations
%                        An integer which describes the number of random
%                        permutations to be used to calculate significance.
%                        Defaults to 10,000.
%                userOptions.distanceMeasure
%                        A string descriptive of the distance measure to be used
%                        to compare two RDMs. Defaults to 'Spearman'.
%
% The following files are saved by this function:
%        userOptions.rootPath/Statistics/
%                userOptions.analysisName_Significance.csv
%                        A .csv (comma-separated value) file with the following
%                        (headed) columns:
%                                First_RDM_name
%                                Second_RDM_name
%                                r
%                                p_perm
%                                p_conv(unreliable?)
%                        and one further row for each comparison in the complete
%                        bipartite graph on all the referenceRDMs and all the
%                        comparedRDMs.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_testSignificance_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%  
% Cai Wingfield 11-2009, 6-2010

function testSignificance(referenceRDMsCell, comparedRDMsCell, userOptions)

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
if ~isfield(userOptions, 'analysisName'), error('testSignificance:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('testSignificance:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'significanceTestPermutations', 10000);
userOptions = setIfUnset(userOptions, 'distanceMeasure', 'Spearman');


StatisticsFileName = [userOptions.analysisName, '_Significance.csv'];
DetailsFileName = [userOptions.analysisName, '_testSignificance_Details.mat'];

% Options for the prompt
promptOptions.functionCaller = 'testSignificance';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Statistics', StatisticsFileName);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFileName);

% Do the prompt
overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag % If files may be (over)written:

	nReferenceRDMStructs = numel(referenceRDMsCell);
	nComparedRDMStructs = numel(comparedRDMsCell);
	
	%% De-cell and concatenate both groups of RDMs
	
	% For reference ones
	for RDMStructI = 1:nReferenceRDMStructs
		thisRDMStruct = referenceRDMsCell{RDMStructI};
		if RDMStructI == 1
			referenceRDMs = thisRDMStruct;
		else
			referenceRDMs = concatenateRDMs(referenceRDMs, thisRDMStruct);
		end%if
	end%for:RDMStructI
	
	% For compared ones
	for RDMStructI = 1:nComparedRDMStructs
		thisRDMStruct = comparedRDMsCell{RDMStructI};
		if RDMStructI == 1
			comparedRDMs = thisRDMStruct;
		else
			comparedRDMs = concatenateRDMs(comparedRDMs, thisRDMStruct);
		end%if
	end%for:RDMStructI

	nPerms = userOptions.significanceTestPermutations;
	nPermsString = num2str(nPerms);

	fprintf(['Testing significance of RDM fits (' nPermsString ' permutations per comparison!)...\n']);

	% Create/open the .csv file containing the statistics
	fileID = fopen(StatisticsFileName, 'wt');

	% Write the column titles to the file
	fprintf(fileID, '%s%s%s%s%s%s%s%s%s', ['First_RDM_name' ',' 'Second_RDM_name' ',' 'r' ',' 'p_perm' ',' 'p_conv(unreliable?)']);
	fprintf(fileID, '\n');
	
	nReferenceRDMs = numel(referenceRDMs);
	nComparedRDMs = numel(comparedRDMs);

	nConditions = size(referenceRDMs(1,1).RDM, 1);

	if nConditions < 8
		fprintf(['(!) There are ' num2str(nConditions) ' (fewer that 8) experimental conditions in this analysis and\n' ...
			'thus using more than ' num2str(nConditions) '! = ' num2str(factorial(nConditions)) ' random permutations will produce unreliable\n' ...
			'results. Instead, each of the ' num2str(factorial(nConditions)) ' permutations will be exhaustively used for best results.' ...
			]);
	end%if

	for referenceRDMi = 1:nReferenceRDMs % For each reference RDMs...
		fprintf(['\t' referenceRDMs(referenceRDMi).name]);
		tic % [Start measuring the time]
		for comparedRDMi = 1:nComparedRDMs % ...and each RoI...
			[r, p, p_conv] = testRDMrelatedness_randomization(referenceRDMs(referenceRDMi).RDM, comparedRDMs(comparedRDMi).RDM, struct('nSignificanceTestPermutations', nPerms, 'corrType', userOptions.distanceMeasure)); % ...compute the pairwise statistics...
			fprintf(fileID, '%s%s%s%s%s%s%s%s%s', [referenceRDMs(referenceRDMi).name ',' comparedRDMs(comparedRDMi).name ',' num2str(r) ',' num2str(p) ',' num2str(p_conv)]); % ...and write the results to the stats file.
			fprintf(fileID, '\n');
			clear p r p_conv;
			fprintf('.');
		end%for
		t = toc; % [stop measuring the time]
		% Verbose outline to the command window:
		fprintf([':\n\t\t' num2str(referenceRDMi) ' of ' num2str(nReferenceRDMs) ': Done.\n']);
		fprintf(['\t\tTime to compute: ~' num2str(ceil(t)) 's.\n']);
		fprintf(['\t\tEstimated time remaining: ~' num2str(ceil((nReferenceRDMs-referenceRDMi)*t)) 's.\n']);
		clear t;
	end%for

	timeStamp = datestr(now);

	fprintf(['Saving statistics data to ' fullfile(pwd, StatisticsFileName) '...\n']);
	fclose(fileID); % Close the stats file
	
	fprintf(['Saving Details to ' fullfile(pwd, DetailsFileName) '\n']);
	save(fullfile(userOptions.rootPath, 'Details',DetailsFileName), 'timeStamp', 'userOptions');

end%if

cd(returnHere);
