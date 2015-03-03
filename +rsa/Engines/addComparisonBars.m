function y = addComparisonBars(pairwisePs,cYMax,threshold)

% this function adds the pairwise copmparison bars to a figure; the 'hold on'
% should have been set before executing this function

% pairWisePs: an nxn pairwise comparison matrix. This is supposed to be
% sorted in descending order of average values (large to small heights)
% a horizontal bar would be displayed whenever the (one-sided) comparison 
% reaches a threshold.
% HN, May 2013
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council
nTestRDMs = size(pairwisePs,1);
eachLineHeight = 0.03;
nSignificantComparison = 0;

nPairwiseLines = (nTestRDMs * (nTestRDMs - 1)) / 2;
padding = 2;
% 	nYMax = ((nPairwiseLines + padding) * eachLineHeight + 1) * cYMax;
nYMax = (nSignificantComparison+1) * cYMax;
    
    
pairwiseI = 0;
padding = 2;
for i = 1:nTestRDMs
    for j = i+1:nTestRDMs
        
        % The p value for i (the lower bar) being significantly smaller than j (the higher bar) is the p value for j being significantly higher than i. This is pairwisePs(j,i) (the lt?)
        thisPairwiseP = pairwisePs(i,j);
        
        if thisPairwiseP <= threshold
            y = cYMax;
            pairwiseI = pairwiseI + 1;
            xx = [i j];
            
            y = y + (cYMax * eachLineHeight * padding / 2); %padding
            y = y + (pairwiseI * eachLineHeight * cYMax); %line separation
            yy = [y y];
            line(xx, yy, 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none', 'Color', [0 0 0]);
            nSignificantComparison = nSignificantComparison + 1;
            
        end%if
    end%for:j
end%for:i