% simulateSPMFiles
%
% [betaCorrespondence_true betaCorrespondence_noisy]...
%                     = simulateSPMFiles(userOptions, simulationOptions)
%
% Creates and saves SPM files for each subject specified in userOptions
% containing data clustered according to simulationOptions, both 'true' and
% 'noisy'.
%
% Cai Wingfield 6-2010

function [betaCorrespondence_true, betaCorrespondence_noisy] = simulateSPMFiles(userOptions, simulationOptions)

returnHere = pwd; % We'll return to the pwd when the function has finished

DetailsFilename = [userOptions.analysisName, '_simulateSPMFiles_Details.mat'];

promptOptions.functionCaller = 'simulateSPMFiles';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag

	nSubjects = numel(userOptions.subjectNames);
	
	% For each subject...
	for subject = 1:nSubjects
	
		thisSubject = userOptions.subjectNames{subject};
	
		fprintf(['Simulating fMRI data for subject ' num2str(subject) ' of ' num2str(nSubjects) ': ' thisSubject '...']);
	
		tic%1
	
		% Get the directory where these 
		saveDir = replaceWildcards(userOptions.betaPath, '[[subjectName]]', thisSubject, '[[betaIdentifier]]', '');
		try
			cd(saveDir);
		catch
			mkdir(saveDir);
			cd(saveDir);
		end%try
	
		% Bet the [nConditions nVoxels]-sized B matrices
		[B_true, B_noisy] = simulateClusteredfMRIData(simulationOptions);
	
		nConditions = size(B_true, 1);
	
		% For each condition...
		for condition = 1:nConditions
	
			% Make the beta images the right size
			trueBetaImage = zeros(simulationOptions.volumeSize_vox);
			noisyBetaImage = zeros(simulationOptions.volumeSize_vox);
	
			% Pull out the right conditions
			trueBetaImage(:) = B_true(condition,:);
			noisyBetaImage(:) = B_noisy(condition,:);
	
			% Sort out the name of these beta images
			betaCorrespondence_true(1,condition).identifier = ['Beta_' userOptions.conditionLabels{condition} '_true.img'];
			betaCorrespondence_noisy(1,condition).identifier = ['Beta_' userOptions.conditionLabels{condition} '_noisy.img'];
	
			% Make an SPM struct to hold the true image
			V = spawnSPMStruct();
			V.fname = fullfile(saveDir, betaCorrespondence_true(1,condition).identifier);
			V.dim = simulationOptions.volumeSize_vox;
	
			% Write this image
			spm_write_vol(V, trueBetaImage);
	
			% Do the same for the noisy image
			V.fname = fullfile(saveDir, betaCorrespondence_noisy(1,condition).identifier);
	
			% Write the image
			spm_write_vol(V, noisyBetaImage);
	
		end%for:condition
	
		t = toc;%1
	
		fprintf(['[' num2str(ceil(t)) 's]\n']);
	
	end%for:subject

	%% Save relevant info

	timeStamp = datestr(now);

	fprintf(['Saving Details to ' fullfile(pwd, DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions', 'betaCorrespondence_true', 'betaCorrespondence_noisy');
	
else
	
	details = load(fullfile(userOptions.rootPath, 'Details', DetailsFilename));
	betaCorrespondence_true = details.betaCorrespondence_true;
	betaCorrespondence_noisy = details.betaCorrespondence_noisy;

end%if

cd(returnHere); % Go back (probably will never have left)
