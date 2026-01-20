function [EEG, badChans] = happe_detectBadChans_listOnly(EEG, params, pt)
% happe_detectBadChans_listOnly
% Detect bad channels using HAPPE logic but DO NOT remove them.
% Instead, return a list of detected bad channels.
%
% Usage:
%   [EEG, badChans] = happe_detectBadChans_listOnly(EEG, params, pt)
%
% Outputs:
%   EEG      - Original EEG structure (unchanged)
%   badChans - struct containing detected bad channels
%              .indices : channel indices
%              .labels  : channel labels
%
% This function mirrors happe_detectBadChans but prevents automatic
% channel rejection by capturing rejection masks instead of applying them.

fprintf('Detecting bad channels (no automatic removal)...\n');

origEEG = EEG;   % preserve original EEG
badIdx = [];

%% ---------------- Low-density pipeline ----------------
if params.lowDensity
    % Run clean_rawdata but do NOT reject
    tmpEEG = pop_clean_rawdata(EEG, ...
        'FlatlineCriterion', 10, ...
        'ChannelCriterion', .5, ...
        'LineNoiseCriterion', 'off', ...
        'Highpass', 'off', ...
        'BurstCriterion', 'off', ...
        'WindowCriterion', 'off', ...
        'BurstRejection', 'off', ...
        'Distance', 'Euclidian', ...
        'RejectChannels', 'off');

    % Channels clean_rawdata WOULD have removed
    if isfield(tmpEEG.etc,'clean_channel_mask')
        badIdx = find(~tmpEEG.etc.clean_channel_mask);
    end

%% ---------------- High-density (HAPPE v2) pipeline ----------------
else
    switch pt
        case 1
            upThresh = 1.8935; chanCrit = 0.485; lnCrit = 7.1;
        case 2
            upThresh = 2.1316; chanCrit = 'off'; lnCrit = 'off';
    end

    % ---- pop_rejchan (spectral) ----
    [~, indelec] = pop_rejchan(EEG, ...
        'elec', 1:EEG.nbchan, ...
        'threshold', [-5 upThresh], ...
        'norm', 'on', ...
        'measure', 'spec', ...
        'freqrange', [1 100], ...
        'reject', 0);   % do NOT remove

    badIdx = unique([badIdx indelec]);

    % ---- clean_rawdata (correlation / line noise) ----
    tmpEEG = pop_clean_rawdata(EEG, ...
        'FlatlineCriterion', 'off', ...
        'ChannelCriterion', chanCrit, ...
        'LineNoiseCriterion', lnCrit, ...
        'Highpass', 'off', ...
        'BurstCriterion', 'off', ...
        'WindowCriterion', 'off', ...
        'BurstRejection', 'off', ...
        'Distance', 'Euclidian', ...
        'RejectChannels', 'off');

    if isfield(tmpEEG.etc,'clean_channel_mask')
        badIdx = unique([badIdx find(~tmpEEG.etc.clean_channel_mask)]);
    end
end

%% ---------------- Output ----------------
badChans.indices = sort(badIdx);
badChans.labels  = {EEG.chanlocs(badIdx).labels};

EEG = origEEG;   % ensure EEG unchanged

fprintf('Detected %d bad channels:\n', numel(badIdx));
if isempty(badIdx)
    fprintf('  NONE\n');
else
    fprintf('  %s\n', strjoin(badChans.labels, ', '));
end

end
