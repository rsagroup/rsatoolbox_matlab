function [varargout] = simulateClusteredfMRIData_fullBrain(simulationOptions)
%
% this function simulates fMRI data with a categorical effect inserted at a
% specified location. The remaining regions contain noise that is temporally
% and spatially smooth yet lacks any particular categorical structure.
% the size of the simulated volume is indicated by the volumeSize_vox field
% and the effect center by effectCen. The volume of the simulated effect
% (within the brain volume) is controlled by volumeSize_vox.
% Hamed Nili 2012
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

nNoisyPatterns = nargout - 2;

%% Generate B patterns
b = generateBetaPatterns(simulationOptions.clusterSpec, prod(simulationOptions.volumeSize_vox));

nVox_wholeBrain = prod(simulationOptions.brainVol);
msk = zeros(simulationOptions.brainVol);
effectCen = simulationOptions.effectCen;
msk([effectCen(1)+1:effectCen(1)+simulationOptions.volumeSize_vox(1)],[effectCen(2)+1:effectCen(2)+simulationOptions.volumeSize_vox(2)],[effectCen(3)+1:effectCen(3)+simulationOptions.volumeSize_vox(3)]) = 1;
B_true = zeros(simulationOptions.nConditions,nVox_wholeBrain);
B_true(:,find(msk)) = b;
nConditions = size(B_true,1);
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

	[E, smoothedYfilename_ignore] = spatiallySmooth4DfMRI_mm(E, simulationOptions.brainVol, simulationOptions.spatiotemporalSmoothingFWHM_mm_s(1:3), simulationOptions.voxelSize_mm);

	% Smooth across time
	E = temporallySmoothTimeSpaceMatrix(E, simulationOptions.spatiotemporalSmoothingFWHM_mm_s(4) / simulationOptions.TR);

	%% Do GLM for Y_true matrix
	Y_true = X * B_true;
    varargout{2} = msk;
    varargout{3} = Y_true;
	%% Do GLM for Y_noise matrix
	Y_noisy = Y_true + E;
	B_noisy = inv(X' * X) * X' * Y_noisy;
    
    fMRI.B = B_noisy;
    fMRI.Y = Y_noisy;
    fMRI.X = X;
    fMRI.groundTruth = b;
       
	varargout{o + 3} = fMRI;
    
	
	clear E Y_noisy B_noisy;
	
end%for:o

end%function
