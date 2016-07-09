function [cognitiveX,BVstimProt,standardIndexSequence,hirf_ms]=generateCognitiveModel_fastButTrialsNeedToStartOnVols(sequence,stimDur_ms,trialDur_ms,nTRvols_final,TRvol_ms,nSkippedVols,monitor,scaleTrialResponseTo1)

% this function generates a design matrix of hemodynamic response
% predictors from a condition sequence.
%
% this function requires the trial duration to be a multiple of the TRvol
% and the trials to start in register with the volume acquisitions.
% cf. related function:
% generateCognitiveModel_generalButSlow.m
%
% NOTES
% highest trial index is assumed to be null (rest or fixation) trial
%
% design matrix is generated in milliseconds and then subsampled at the TRvol
% resolution
%
% EDIT: CW 6-2010: No longer saves 'cognitiveX.txt'
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

%% control variables
if ~exist('scaleTrialResponseTo1','var');
    scaleTrialResponseTo1=1;
end

% if ~scaleTrialResponseTo1
%     button=questdlg('Do you want to scale the trial response to peak at 1 (in ms resolution)?')
%     if strcmp(button,'Yes')
%         scaleTrialResponseTo1=1;
%     elseif strcmp(button,'Cancel')
%         error('design generation aborted by user')
%     end
% end

if ~exist('monitor','var')
    monitor=1;
end

trialDur_TRvol=trialDur_ms/TRvol_ms;
if trialDur_TRvol~=round(trialDur_TRvol)
    error('ERROR: trialDur_ms must be a multiple of TRvol_ms.');
end

% make sequence a column
sequence=sequence(:);


%% give trial types consecutive integer values
% existingTrialTypeIDsVector=zeros(max(sequence),1);
% existingTrialTypeIDsVector(sequence)=1;
% existingTrialTypeIDs=find(existingTrialTypeIDsVector);
existingTrialTypeIDs=unique(sequence);

standardIndexSequence=zeros(size(sequence));
nTrialTypes=length(existingTrialTypeIDs); %including null trials (highest trial-type index)
for trialTypeI=1:nTrialTypes
    standardIndexSequence(find(sequence==existingTrialTypeIDs(trialTypeI)))=trialTypeI;
end

[trialTypeIForTrialI,index]=sort(standardIndexSequence);

nTrials=length(standardIndexSequence);
%nTRvolsPerTrial=trialDur_ms/TRvol_ms;

i=1:nTrials;
startStim(i)=(i-1)*trialDur_ms;
stopStim=startStim+stimDur_ms;



%% determine start and end times of stimuli in ms and display BV-style stimprot

BVstimProt=[startStim(index); stopStim(index)]';
BVstimProt=BVstimProt-nSkippedVols*TRvol_ms;

if monitor
    disp('_____________________________________________________________________');
    disp('trials'' start and end times in ms (as used in brainvoyager prt files');
    BVstimProt
    disp('_____________________________________________________________________');
end


%% determine hemodynamic response at TR resolution
rectangularNeuralResponse_ms=ones(1,stimDur_ms);
hirf_ms=boyntonModel(1);
hirf_ms=hirf_ms/max(hirf_ms);

hemodynamicTrialResponse_ms=conv(rectangularNeuralResponse_ms,hirf_ms);

if scaleTrialResponseTo1
    hemodynamicTrialResponse_ms=hemodynamicTrialResponse_ms/max(hemodynamicTrialResponse_ms);
    % scale to have a maximum of 1 (at ms resolution, final max at TR
    % resolution (likely lower than 1) will reflect, by how much the acquired
    % vols miss the max HR)
end

lengthHR_ms=20*1000; % 20 seconds
lengthHR_TRvol=lengthHR_ms/TRvol_ms;

%samplingTimes_ms=[0:lengthHR_TRvol-1].*TRvol_ms+TRvol_ms/2; % sample in the middle of each volume acquisition
samplingTimes_ms=[0:lengthHR_TRvol-1].*TRvol_ms; % sample at the beginning of each volume acquisition (identical to BV)

hemodynamicTrialResponse_TRvol=hemodynamicTrialResponse_ms(samplingTimes_ms+1);

if monitor
    pageFigure(100); clf;
    
    subplot(2,1,1);
    plot(hirf_ms,'k:'); hold on;
    plot(hemodynamicTrialResponse_ms,'g');
    plot(samplingTimes_ms+1,hemodynamicTrialResponse_TRvol,'r');
    xlabel('time [ms]');
    title({'\bfBOLD hemodynamic response functions', '\rmimpulse (black dotted), trial-event (green), TR-subsampled trial-event (red)'});

    subplot(2,1,2);
    periodogram(hirf_ms);
    
    
    
    figure(90); hold on;
    h=plot(hemodynamicTrialResponse_ms,'Color',[0 .5 1],'LineWidth',3);
    legend(h,'BOLD');
end


%% make trial-impulse design matrix at TR resolution
durationX_TRvol=nTrials*trialDur_TRvol-nSkippedVols;
nPreds=nTrialTypes-1;
trialImpulseX_TRvol=zeros(durationX_TRvol,nPreds);
for trialI=1:nTrials
    if standardIndexSequence(trialI)~=nTrialTypes % if not a null trial
        if (trialI-1)*trialDur_TRvol+1-nSkippedVols>0
            trialImpulseX_TRvol((trialI-1)*trialDur_TRvol+1-nSkippedVols,standardIndexSequence(trialI))=1;
        else
            warning('WARNING: nSkippedVols specified to skip non-null trials.');
        end
    end
end

if monitor
    figure(200); clf;
    plot(trialImpulseX_TRvol);
    title('trial-impulse design matrix at TR resolution');
end


%% convolve with the boynton-model of the hemodynamic response at TR resolution

for predI=1:nPreds
    cognitiveX(:,predI)=conv(trialImpulseX_TRvol(:,predI),hemodynamicTrialResponse_TRvol); % c)urrent pred)ictor of h)emodynamic r)esponse (in TRvols)
end

if ~isempty(nTRvols_final)
    cognitiveX=cognitiveX(1:nTRvols_final,:);
end

if monitor
    figw(300); clf;
    plot(cognitiveX,'LineWidth',2);
    xlabel(['time [TRvol=',num2str(TRvol_ms),'ms]']);
    title('\bfBOLD: design matrix of hemodynamic response predictors');
end

%save('cognitiveX.txt','cognitiveX','-ASCII')

if monitor==2
    evaluateDesign(cognitiveX)
end

end%function
