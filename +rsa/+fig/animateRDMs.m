function animateRDMs(RDMs, output_filename, delay_time, loop_count, rankTransform01, clims, showColorbar, aspect, colorScheme)
% EXAMPLE USAGE:
%
% for i = 1:30
%    RDMs(i).RDM = pdist(randn(10,10));
%    RDMs(i).name = ['RDM ', num2str(i)];
% end
%
% animateRDMs(RDMs, 'animated.gif', 0.1, Inf, true, [0 1], true, 1, 'jet');
%
% Cai Wingfield 2015-02


    %% define default behavior
    if ~exist('rankTransform01','var'), rankTransform01=true; clims=[0 1]; end
    if ~exist('clims','var'), clims=[]; end
    if ~exist('showColorbar','var'), showColorbar=true; end
    if ~exist('aspect', 'var') || isempty(aspect), aspect = 2/3; end

    if ~exist('delay_time','var'), delay_time = 0; end
    if ~exist('loop_count','var'), loop_count = Inf; end

    %% display RDMs

    % get current highest figure number so we don't replace anything.
    figHandles = findall(0, 'Type', 'figure');
    maxFigHandle = 0;
    for fig_i = 1 : length(figHandles)
        maxFigHandle = max(maxFigHandle, figHandles(fig_i).Number);
    end%for:fig_i
    newFigHandle = maxFigHandle +  1;

    % show each RDM in turn
    % We assume that RDMs are passed in in a 1-d struct for now.
    image_stack = NaN;

    % Initial value of map
    % Hooray for Matlab's weird typing!
    map = 256;

    for RDM_i = 1 : length(RDMs)
        this_fig_handle = newFigHandle + RDM_i;
        showRDMs(RDMs(RDM_i), this_fig_handle, rankTransform01, clims, showColorbar, aspect, [], colorScheme);
        f = getframe(gcf);
        [im, map] = rgb2ind(f.cdata, map, 'nodither');

        % put it on the stack
        if isnan(image_stack)
            image_stack = im;
        else
            image_stack = cat(4, image_stack, im);
        end%if

        % Close this figure
        close;
    end%for:RDMs

    %% Now save the anigif

    imwrite(image_stack, map, output_filename, 'DelayTime', delay_time, 'LoopCount', loop_count);

end%function
