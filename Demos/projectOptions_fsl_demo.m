function userOptions = projectOptions_fsl_demo()
% projectOptions_fsl_demo is an options file for the FSL demo tutorial
%
% Roni Maimon 9-2017 
%__________________________________________________________________________
% Copyright (C)

%% Project details
% This name identifies a collection of files which all belong to the same run of a project.
userOptions.analysisName = 'DEMO_FSL'; % this is renamed in the code for demos 3-4.

% This is the root directory of the project.
userOptions.rootPath = [pwd,filesep,'DEMO_FSL'];


%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL SETUP %%
%%%%%%%%%%%%%%%%%%%%%%%%

% The list of subjects to be included in the study.
userOptions.subjectNames = { ...
	'T05', ...
	'T06', ...
	};

%% FSL Specific parameters
userOptions.run_names = { ...
    'A1b','A2b','A3b','A4b' ...
    };
userOptions.featsPrefix = '5_';
userOptions.featsSuffix = '';
userOptions.featsPath = [pwd,filesep,'FSLData',filesep,'[[subjectName]]',...
                        filesep,'models',filesep,'glm',filesep,'fingers', ...
                        filesep, '[[featPrefix]]','[[runName]]','[[featSufix]]','.feat'];

userOptions.copes = { ...
    1,2,3,4,5 ...
    };

%% End of FSL specific parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS PREFERENCES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First-order analysis
% Text lables which may be attached to the conditions for MDS plots.
userOptions.conditionLabels = { ...
	'F1', ...
	'F2', ...
    'F3', ...
	'F4', ...
    'F5'
	};

userOptions.useAlternativeConditionLabels = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEATUERS OF INTEREST SELECTION OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%% %% %% %% %%
	%% fMRI  %% Use these next three options if you're working in fMRI native space:
	%% %% %% %% %%
	
	% The path to a stereotypical mask data file is stored (not including subject-specific identifiers).
	% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
	% "[[maskName]]" should be used as a placeholder to denote an entry in userOptions.maskNames
	userOptions.maskPath = [pwd,filesep,'FSLData',filesep,'[[subjectName]]',filesep,'models',filesep,'rsa',filesep,'masks',filesep,'[[maskName]].nii.gz'];
		
		% The list of mask filenames (minus .hdr extension) to be used.
	userOptions.maskNames = { ...
			'S1_handknob_alldigits'...
			};
end%function
