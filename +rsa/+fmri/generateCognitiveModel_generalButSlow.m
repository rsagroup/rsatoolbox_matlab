function [cognitiveX,BVstimProt]=generateCognitiveModel_generalButSlow(sequence,stimDur_ms,trialDur_ms,nTRvols_final,TRvol_ms,nSkippedVols,monitor)


import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

% this function generates a design matrix of hemodynamic response
% predictors from a condition sequence.
%
% cf. related function:
% generateCognitiveModel_fastButTrialsNeedToStartOnVols.m
%
% NOTES
% highest trial index is assumed to be null (rest or fixation) trial
%
% design matrix is generated in milliseconds and then subsampled at the TRvol
% resolution

if ~exist('monitor','var')
    monitor=1;
end


% make sequence a column
sequence=sequence(:)


%% determine start and end times of stimuli in ms

% give trial types consecutive integer values
existingTrialTypeIDsVector=zeros(max(sequence),1);
existingTrialTypeIDsVector(sequence)=1;
existingTrialTypeIDs=find(existingTrialTypeIDsVector);
newSequence=zeros(size(sequence));
nTrialTypes=length(existingTrialTypeIDs); %including null trials (highest trial-type index)
for trialTypeI=1:nTrialTypes
    newSequence(find(sequence==existingTrialTypeIDs(trialTypeI)))=trialTypeI;
end
sequence=newSequence;


nTrials=size(sequence,1);
%nTRvolsPerTrial=trialDur_ms/TRvol_ms;

i=1:nTrials;
startStim(i)=(i-1)*trialDur_ms;
stopStim=startStim+stimDur_ms;

[trialTypeIForTrialI,index]=sort(sequence);

BVstimProt=[startStim(index); stopStim(index)]';

BVstimProt=BVstimProt-nSkippedVols*TRvol_ms;


if monitor
    disp('_____________________________________________________________________');
    disp('trials'' start and end times in ms (as used in brainvoyager prt files');
    BVstimProt
    disp('_____________________________________________________________________');
end


%% make rectangular stimulation predictors
durationX_ms=nTrials*trialDur_ms-nSkippedVols*TRvol_ms;
nPreds=nTrialTypes-1;
cognitiveNeuralX_ms=zeros(durationX_ms,nPreds);
for trialI=1:nTrials
    if trialTypeIForTrialI(trialI)~=nTrialTypes
        cognitiveNeuralX_ms(BVstimProt(trialI,1)+1:BVstimProt(trialI,2),trialTypeIForTrialI(trialI))=1;
    end
end



%% convolve with the boynton-model of the hemodynamic response and subsample at TRvol resolution
cognitiveX=zeros(nTRvols_final,nPreds);

hirf_ms=boyntonModel(1);   % Hemodynamic Impuls Response Function (resolution = 1 ms)
hirf_ms=hirf_ms/max(hirf_ms); % scale to have a maximum of 1

% determine single isolated trial response, such that it can be scaled to peak at 1
rectangularNeuralResponse_ms=ones(1,stimDur_ms);
hemodynamicTrialResponse_ms=conv(rectangularNeuralResponse_ms,hirf_ms);

figure(101); clf; hold on;
plot(hirf_ms,':k');
plot(hemodynamicTrialResponse_ms,'g');
title('hemodynamic response functions: impulse (black dotted), trial-event (green)');



%samplingTimes_ms=[0:nTRvols_final-1].*TRvol_ms+TRvol_ms/2; % sample in the middle of each volume acquisition
samplingTimes_ms=[0:nTRvols_final-1].*TRvol_ms; % sample at the beginning of each volume acquisition (identical to BV)

mx_ms=max(samplingTimes_ms)+1;

hp=progressMonitor(1,size(cognitiveNeuralX_ms,2),'generating design matrices at ms resolution...');
for predI=1:size(cognitiveNeuralX_ms,2)
    progressMonitor(predI,size(cognitiveNeuralX_ms,2),'generating design matrices at ms resolution...',hp);
    cHRpred_ms=conv(cognitiveNeuralX_ms(:,predI),hirf_ms/max(hemodynamicTrialResponse_ms)); % c)urrent pred)ictor of h)emodynamic r)esponse (in ms)
    %cognitiveX_ms(:,predI)=cHRpred_ms;
    cHRpred_ms=[cHRpred_ms;zeros(mx_ms-numel(cHRpred_ms),1)];
    cognitiveX(:,predI)=cHRpred_ms(samplingTimes_ms+1);
end
closeProgressMonitor(hp);

save('cognitiveX1.txt','cognitiveX','-ASCII')
save('cognitiveX','cognitiveX')

if monitor
    figure(301); clf; 
    subplot(2,1,1); hold on;
    %plot(cognitiveX_ms, ':'); % in ms
    plot(cognitiveNeuralX_ms,'k');
    plot(samplingTimes_ms+1,cognitiveX); % in TRvols

    subplot(2,1,2); hold on;
    plot(cognitiveX);
end

if monitor==2
    evaluateDesign(cognitiveX)
end

    
    

