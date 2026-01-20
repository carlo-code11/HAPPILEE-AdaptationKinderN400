function plot_erp_allchans(EEG_S1, EEG_S2, condLabels, figTitle)
% plot_erp_allchans_layout
% Plot ERPs for two conditions for all channels arranged by scalp layout
%
% Inputs:
%   EEG_S1, EEG_S2 : epoched EEGLAB structs (same chans & times)
%   condLabels     : {'cond1','cond2'}
%   figTitle       : figure title

    % --- sanity checks ---
    assert(EEG_S1.nbchan == EEG_S2.nbchan, 'Channel count mismatch');
    assert(isequal({EEG_S1.chanlocs.labels}, {EEG_S2.chanlocs.labels}), ...
        'Channel labels do not match');

    % --- compute ERPs ---
    ERP1 = mean(EEG_S1.data, 3);   % chan × time
    ERP2 = mean(EEG_S2.data, 3);

    t = EEG_S1.times;

    % --- get 2D channel coordinates ---
%     radius = [EEG_S1.chanlocs.radius];
%     radius = [radius(1:15) radius(18)];
%     theta = [EEG_S1.chanlocs.theta];
%     theta = [theta(1:15) theta(18)];
    x = [EEG_S1.chanlocs.radius] .* sind([EEG_S1.chanlocs.theta]);
    y = [EEG_S1.chanlocs.radius] .* cosd([EEG_S1.chanlocs.theta]);
    

    

    % normalize to [0,1] figure space
    x = rescale(x, 0.1, 0.9);
    y = rescale(y, 0.1, 0.9);

    % --- figure ---
    figure('Color','w','Units','normalized','Position',[0.1 0.1 0.8 0.8]);
%     sgtitle(figTitle);

    % --- plot each channel ---
    for ch = 1:EEG_S1.nbchan - 2

        ax = axes('Position', [ ...
            x(ch) - 0.03, ...
            y(ch) - 0.03, ...
            0.08, ...
            0.08 ]);

        plot(t, ERP1(ch,:), 'LineWidth', 1.2); hold on;
        plot(t, ERP2(ch,:), 'LineWidth', 1.2);

        set(ax, ...
            'YDir', 'reverse', ...
            'Box', 'off');

        xlim([t(1) t(end)]);
        ylim([-5 5]);
        if ch < EEG_S1.nbchan - 2
            title(EEG_S1.chanlocs(ch).labels, ...
                  'FontSize', 7);
        else
            title('CZ', ...
                  'FontSize', 7);
        end
    end

    % --- legend ---
    legend(condLabels, ...
           'Orientation','horizontal', ...
           'Position',[0.1 0.01 0.15 0.03]);

end
