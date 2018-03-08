% this function adds the pairwise copmparison bars to a figure; the 'hold on'
% should have been set before executing this function

% pairWisePs: an nxn pairwise comparison matrix. This is supposed to be
% sorted in descending order of average values (large to small heights)
% a horizontal bar would be displayed whenever the (one-sided) comparison 
% reaches a threshold.
% HN, May 2013
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

function y = addComparisonBars(pairwisePs,cYMax,threshold)

    import rsa.*
    import rsa.fig.*
    import rsa.rdm.*
    import rsa.stat.*
    import rsa.util.*


    nTestRDMs = size(pairwisePs,1);
    eachLineHeight = 0.03;

    padding = 2;

    % Set `y` here in case no comparisons are significant.
    y = cYMax;


    nSignificantComparison = 0;
    pairwiseI = 0;
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
end%function
