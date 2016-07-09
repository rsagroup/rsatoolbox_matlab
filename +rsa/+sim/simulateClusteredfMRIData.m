function [varargout] = simulateClusteredfMRIData(simulationOptions)
% simulateClusteredfMRIData
%
% [B_true, B_noisy[, B_noisy2[, ...]]] = simulateClusteredfMRIData(simulationOptions)
%
% Simulate both a 'true' beta matrix with clustered data (according to
% simulationOptions) and a 'noisy' one contaminated with simulated scanner
% noise.
%
% Based on code by Niko and Hamed
%
% Cai Wingfield 6-2010, 8-2010
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

nNoisyPatterns = nargout - 1;

%% Generate B patterns
B_true = generateBetaPatterns(simulationOptions.clusterSpec, prod(simulationOptions.volumeSize_vox));

varargout{1} = B_true;

nConditions = size(B_true, 1);
nVoxels = size(B_true,2);

%% Generate X
sequence = [repmat([1:nConditions],1,simulationOptions.nRepititions),(nConditions+1)*ones(1,ceil(nConditions*simulationOptions.nRepititions/3))];
sequence = randomlyPermute(sequence);
nTrials = numel(sequence);
nTRvols = (simulationOptions.trialDuration/simulationOptions.TR)*nTrials;
nSkippedVols = 0;
monitor = 0;
scaleTrialResponseTo1 = 1;

[X,BV_ignore,standardIndexSequence_ignore,hirf_ms_ignore] = generateCognitiveModel_fastButTrialsNeedToStartOnVols(sequence,simulationOptions.stimulusDuration*1000,simulationOptions.trialDuration*1000,nTRvols,simulationOptions.TR*1000,nSkippedVols,monitor,scaleTrialResponseTo1);

nTimePoints = size(X,1);
sig = sqrt(simulationOptions.scannerNoiseLevel);

varargout{1} = B_true;

for o = 1:nNoisyPatterns
	
	%% Generate E matrix
	E = randn(nTimePoints, nVoxels);
	E = sig * E;

	[E, smoothedYfilename_ignore] = spatiallySmooth4DfMRI_mm(E, simulationOptions.volumeSize_vox, simulationOptions.spatiotemporalSmoothingFWHM_mm_s(1:3), simulationOptions.voxelSize_mm);

	% Smooth across time
	E = temporallySmoothTimeSpaceMatrix(E, simulationOptions.spatiotemporalSmoothingFWHM_mm_s(4) / simulationOptions.TR);

	%% Do GLM for Y_true matrix
	Y_true = X * B_true;

	%% Do GLM for Y_noise matrix
	Y_noisy = Y_true + E;
	B_noisy = inv(X' * X) * X' * Y_noisy;
    
    fMRI.B = B_noisy;
    fMRI.Y = Y_noisy;
    fMRI.X = X;
       
	varargout{o + 1} = fMRI;
	
	clear E Y_noisy B_noisy;
	
end%for:o

end%function
