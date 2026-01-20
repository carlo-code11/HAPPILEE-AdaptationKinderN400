% setParams() - A helper function for HAPPE that prints
%                out the user-entered parameters to the command window in a
%                way that is easy to read. Allows the user to check their
%                parameters prior to confirming them.
%
% Usage: 
%   >> params = setParams(params, pre_exist, reprocessing, changed_params, happeDir)
%
% Inputs:
%   params         - A struct containing all of the parameters needed to 
%                    run HAPPE.
%   pre_exist      - A flag indicating whether or not these are pre-existing
%                    paramaters (0 for no, 1 for yes).
%   reprocessing   - A flag indicating whether or not the data is being
%                    reprocessed (0 for no, 1 for yes).
%   changed_params - A flag indicating whether or not the parameters are
%                    being changed.
%   happeDir       - The directory containing all the HAPPE scripts.
%                    Included for testing purposes.
%
% Outputs:
%   params         - A struct containing all the parameters needed to run
%                    HAPPE
%
% NOTE: Changes to the way any parameters are encoded, or changes in the
% number of parameters (adding/subtracting) will result in this code
% needing to be updated to reflect those changes. There is no easy way to
% make it reliant on other scripts... Sorry :(
%
% Author: A.D. Monachino, PINE Lab at Northeastern University, 2022
%
% This file is part of HAPPE.
% Copyright 2018, 2021 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% HAPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% HAPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with HAPPE. If not, see <https://www.gnu.org/licenses/>.

function params = setParams(params, preExist, reprocessing, changedParams, happeDir)
paramChoice = 'na' ;
while true
    %% BREAK IF NOT CHANGING ANY PRE-EXISTING PARAMETERS
    if preExist && ~changedParams; break; end
    %% IF CHANGING PARAMETERS
    % SET THE CHANGE MESSAGE BASED ON THE AVAILABLE OPTIONS ACCORDING TO
    % THE PARADIGM AND REPROCESSING STATE
    changeMessage = [''] ;
    if reprocessing
        if isfield(params, 'paradigm')
            if params.paradigm.ERP.on
                changeMessage = ['filter, segmentation, interpolation, ' ...
                    'segment rejection, re-referencing, save format, ' ...
                    'visualizations' changeMessage] ;
            else
                changeMessage = ['muscIL, segmentation, interpolation, ' ...
                    'segment rejection, re-referencing, save format, ' ...
                    'visualizations' changeMessage] ;
            end
        end
    else
        if isfield(params, 'paradigm')
            if params.paradigm.task
                if params.paradigm.ERP.on
                    changeMessage = ['onset tags, file import information, ' ...
                        'channels of interest, channel renaming, line noise ' ...
                        'frequency, line noise reduction, resampling, ' ...
                        'filter, bad channel detection, wavelet thresholding, ' ...
                        'segmentation, '...
                        'interpolation, segment rejection, ' ...
                        're-referencing, save format, visualizations' ...
                        changeMessage] ;
                else
                    changeMessage = ['onset tags, file import information, ' ...
                        'channels of interest, channel renaming, line noise ' ...
                        'frequency, line noise reduction, resampling, ' ...
                        'filter, bad channel detection, wavelet thresholding, ' ...
                        'muscIL, segmentation, ' ...
                        'interpolation, segment rejection, ' ...
                        're-referencing, save format, visualizations' ...
                        changeMessage] ;
                end
            else
                changeMessage = ['file import information, channels of ' ...
                    'interest, line noise frequency, line noise reduction,' ...
                    ' resampling, filter, bad channel detection, wavelet ' ...
                    'thresholding, muscIL, ' ...
                    'segmentation, interpolation segment rejection, ' ...
                    're-referencing, save format, visualizations' ...
                    changeMessage] ;
            end
        end
    end

    %% LIST CHANGEABLE PARAMETERS:
    if changedParams
%         fprintf(['Parameter to change: ' changeMessage]) ;
%         paramChoice = input('> ', 's') ;
        paramChoice = inputdlg(['Parameter to change: ' changeMessage]); 
    end
    
    if ~preExist
        %% DENSITY:
        % Set the density as either low or high using user input via the
        % command line. Low density is defined as 30 or fewer channels for
        % HAPPE.
%         fprintf(['Low density data? [Y/N]\nFor HAPPE, low density data ' ...
%             'contains 1 to ~32 channels.\n']) ;
%         params.lowDensity = choose2('N', 'Y') ;
        params.lowDensity = questdlg('Low density data? (For HAPPE, low density data contains 1 to ~32 channels)', ...
                        'Choose Option', ...
                        'yes', 'no', 'yes');
        switch params.lowDensity
            case 'no'
                    params.lowDensity = 0;
            case 'yes'
                    params.lowDensity = 1;
        end                    


        %% PARADIGM (REST VS TASK & ERP):
        % DETERMINE WHETHER THE PARADIGM IS REST OR TASK-RELATED: use user
        % input from the command window. If an invalid input is entered,
        % re-issue prompt until a valid input is entered.
        params.paradigm = struct() ;
%         fprintf(['Data type:\n  rest = Resting-state/baseline EEG\n  task = ' ...
%             'Task-related EEG\n']) ;
%         params.paradigm.task = choose2('rest', 'task') ;
        params.paradigm.task = questdlg('Data type: rest = Resting-state/baseline EEG, task = Task-related EEG', ...
                        'Choose Option', ...
                        'task', 'rest', 'task');
        switch params.paradigm.task
            case 'rest'
                    params.paradigm.task = 0;
            case 'task'
                    params.paradigm.task = 1;
        end  
        % IF TASK-RELATED DATA, PROMPT IF ERP ANALYSIS: use user input from
        % the command window to determine if pre-processing ERP data. If
        % invalid input, repeat until valid input is entered.
        if params.paradigm.task
%             fprintf(['Performing event-related potential (ERP) analysis? ' ...
%                 '[Y/N]\n']) ;
%             params.paradigm.ERP.on = choose2('n', 'y') ;
            params.paradigm.ERP.on = questdlg('Performing event-related potential (ERP) analysis?', ...
                        'Choose Option', ...
                        'yes', 'no', 'yes');
            switch params.paradigm.ERP.on
                case 'no'
                        params.paradigm.ERP.on = 0;
                case 'yes'
                        params.paradigm.ERP.on = 1;
            end  
        else; params.paradigm.ERP.on = 0 ;
        end
        % SET QC FREQS BASED ON THE PARADIGM: Use the paradigm to determine
        % which set of frequencies to use in evaluating pipeline metrics.
        if params.paradigm.ERP.on; params.QCfreqs = [.5, 1, 2, 5, 8, 12, ...
                20, 30, 45, 70] ;
        else; params.QCfreqs = [2, 5, 8, 12, 20, 30, 45, 70] ;
        end
    end
    
    %% TASK ONSET TAGS
%     if params.paradigm.task && (~preExist || (strcmpi(paramChoice, ...
%             'onset tags') && ~reprocessing))
%         fprintf(['Enter the task onset tags, one at a time, pressing ' ...
%                 'enter/return between each entry.\nWhen you have entered ' ...
%                 'all tags, input "done" (without quotations).\n']) ;
%         params.paradigm.onsetTags = UI_cellArray(1,{}) ;
%         if size(params.paradigm.onsetTags,2) > 1
%             fprintf(['Do multiple onset tags belong to a single condition?' ...
%                 ' [Y/N]\nExample: "happy_face" and "sad_face" belong to ' ...
%                 '"faces".\n']) ;
%             params.paradigm.conds.on = choose2('N', 'Y') ; 
%             if params.paradigm.conds.on
%                 params.paradigm.conds.groups = [] ;
%                 fprintf(['Enter the conditions and included onset tags:\n' ...
%                 'Enter each condition as a list, each element seperated by ' ...
%                 'a blank space,\nwith the condition name as the first item ' ...
%                 'in the list. Press enter/return between entries.\nWhen you ' ...
%                 'have entered all conditions, input "done" (without quotations).\n' ...
%                 'Example: faces happy_face sad_face angry_face\n']) ;
%                 while true
%                     temp = split(input('> ', 's'))' ;
%                     if size(temp,2) == 1 && strcmpi(temp, 'done'); break ;
%                     elseif size(temp,2) > 2
%                         diff = size(temp,2) - size(params.paradigm.conds.groups,2) ;
%                         if diff > 0 && size(params.paradigm.conds.groups, 2) ~= 0
%                             params.paradigm.conds.groups = [params.paradigm.conds.groups ...
%                                 cell(size(params.paradigm.conds.groups, ...
%                                 1), diff)] ;
%                         elseif diff < 0 && size(params.paradigm.conds.groups, 2) ~= 0
%                             temp = [temp cell(1, abs(diff))] ;                %#ok<*AGROW>
%                         end
%                          params.paradigm.conds.groups = ...
%                                 [params.paradigm.conds.groups; temp] ;
%                     else
%                         fprintf(['Invalid input: enter a condition name ' ...
%                             'and at least two entries or "done" ' ...
%                             '(without quotations).\n']) ;
%                         continue ;
%                     end
%                 end
%             else; params.paradigm.conds.on = 0 ;
%             end
%         else; params.paradigm.conds.on = 0 ;
%         end
%     end
    
if params.paradigm.task && (~preExist || (strcmpi(paramChoice, 'onset tags') && ~reprocessing))
    
    % Prompt for onset tags using input dialog repeatedly until "done"
    params.paradigm.onsetTags = {};
    while true
        prompt = 'Enter a task onset tag (or type "done" to finish):';
        dlgTitle = 'Task Onset Tags';
        answer = inputdlg(prompt, dlgTitle, [1 50]);
        if isempty(answer)
            % User cancelled
            break;
        end
        tag = strtrim(answer{1});
        if strcmpi(tag, 'done')
            break;
        elseif ~isempty(tag)
            params.paradigm.onsetTags{end+1} = tag; %#ok<*AGROW>
        end
    end

    % If multiple tags exist, ask if they belong to a single condition
    if numel(params.paradigm.onsetTags) > 1
        choice = questdlg(['Do multiple onset tags belong to a single condition?' ...
            '\nExample: "happy_face" and "sad_face" belong to "faces".'], ...
            'Multiple Onset Tags', 'Yes', 'No', 'No');
        params.paradigm.conds.on = strcmpi(choice, 'Yes');

        if params.paradigm.conds.on
            params.paradigm.conds.groups = {};
            % Prompt for condition groups
            while true
                prompt = ['Enter a condition and its onset tags as a space-separated list, ' ...
                          'with the condition name first (e.g., "faces happy_face sad_face"). ' ...
                          'Click Cancel or leave empty when done:'];
                dlgTitle = 'Enter Condition Group';
                answer = inputdlg(prompt, dlgTitle, [1 100]);
                if isempty(answer) || all(cellfun(@isempty, answer))
                    break; % done
                end
                temp = strsplit(strtrim(answer{1}));
                if numel(temp) < 2
                    uiwait(msgbox('Invalid input: enter a condition name and at least one onset tag.','Error','error'));
                    continue;
                end
                params.paradigm.conds.groups = [params.paradigm.conds.groups; temp]; %ok<AGROW>
            end
        else
            params.paradigm.conds.on = 0;
        end
    else
        params.paradigm.conds.on = 0;
    end
end


     %% IMPORTING FILE INFORMATION
    if ~preExist || (strcmpi(paramChoice, 'file import information') && ~reprocessing)
        params.loadInfo = struct() ;
        params.loadInfo = loadingInfo(params, happeDir) ;
    end
    
    %% CHANNEL RENAMING
    if ~preExist || (strcmpi(paramChoice, 'channel renaming') && ~reprocessing)
       params.channelrenameN400.on = questdlg('Do you want to rename channels (to V-,V+,H-,H+,A1,A2)?', ...
                    'Choose Option', ...
                    'yes', 'no', 'yes');
        switch params.channelrenameN400.on
            case 'no'
                    params.channelrenameN400.on = 0;
            case 'yes'
                    params.channelrenameN400.on = 1;
        end 
    end
    
    %% CHANNELS OF INTEREST
    if ~preExist || (strcmpi(paramChoice, 'channels of interest') && ~reprocessing)
        if ~params.loadInfo.chanlocs.inc
            params.chans.IDs = {} ; params.chans.subset = 'all' ;
        else; [params.chans.subset, params.chans.IDs] = determ_chanIDs() ;
        end
    end
    
    %% LINE NOISE FREQUENCY AND METHOD
    if (~preExist || strcmpi(paramChoice, 'line noise frequency')) && ~reprocessing
%         params.lineNoise.freq = input(['Frequency of electrical (line) ' ...
%              'noise in Hz:\nUSA data probably = 60; Otherwise, probably = ' ...
%              '50\n> ']) ; 
        params.lineNoise.freq = inputdlg('Frequency of electrical (line) noise in Hz: USA data probably = 60; Otherwise, probably = 50');
        params.lineNoise.freq = str2num(params.lineNoise.freq{1});

        params.lineNoise.neighbors = [params.lineNoise.freq-10, params.lineNoise.freq-5, ...
            params.lineNoise.freq-2, params.lineNoise.freq-1, params.lineNoise.freq, ...
            params.lineNoise.freq+1, params.lineNoise.freq+2, params.lineNoise.freq+5, ...
            params.lineNoise.freq+10] ;
        
%         fprintf(['Are there any additional frequencies, (e.g., harmonics) to' ...
%             ' reduce? [Y/N]\n']) ;
%         params.lineNoise.harms.on = choose2('n', 'y') ;
        params.lineNoise.harms.on = questdlg('Are there any additional frequencies, (e.g., harmonics) to reduce?', ...
                    'Choose Option', ...
                    'yes', 'no', 'yes');
        switch params.lineNoise.harms.on
            case 'no'
                    params.lineNoise.harms.on = 0;
            case 'yes'
                    params.lineNoise.harms.on = 1;
        end 
        if params.lineNoise.harms.on
%             fprintf(['Enter the frequencies, one at a time, to pass through' ...
%                 ' noise reduction.\nWhen finished entering frequencies, enter ' ...
%                 '"done" (without quotations).\n']) ;
%             indx = 1 ;
%             while true
%                 user_input = input('> ', 's') ;
%                 if strcmpi(user_input, 'done')
%                     params.lineNoise.harms.freqs = unique(params.lineNoise.harms.freqs, ...
%                         'stable') ;
%                     break ;
%                 else
%                     params.lineNoise.harms.freqs(indx) = str2double(user_input) ;
%                     indx = indx + 1 ;
%                 end
%             end
            freqs = inputdlg('Enter the frequencies, seperated by a blank space, to pass through noise reduction:') ;
            params.lineNoise.harms.freqs = str2num(freqs{1});
        else
            params.lineNoise.harms.freqs = [] ;
        end
    end

    if (~preExist || strcmpi(paramChoice, 'line noise reduction')) && ~reprocessing
%         fprintf(['Line noise reduction method:\n  cleanline = Use Tim ' ...
%             'Mullen''s CleanLine\n  notch = Use a notch filter.\n']) ;
%         params.lineNoise.cl = choose2('notch', 'cleanline') ;
        params.lineNoise.cl = questdlg('Line noise reduction method: cleanline = Use Tim Mullens CleanLine, notch = Use a notch filter', ...
                    'Choose Option', ...
                    'cleanline', 'notch', 'cleanline');
        switch params.lineNoise.cl
            case 'notch'
                    params.lineNoise.cl = 0;
            case 'cleanline'
                    params.lineNoise.cl = 1;
        end 
        if params.lineNoise.cl
%             fprintf(['Use legacy or default line noise reduction?\n  default - Default method' ...
%                 ' optimized in HAPPE v2\n  legacy - Method from HAPPE v1 (NOT' ...
%                 ' RECOMMENDED\n']) ;
%             params.lineNoise.legacy = choose2('default', 'legacy') ;
           params.lineNoise.cl = questdlg('Use legacy or default line noise reduction? default - Default method optimized in HAPPE v2, legacy - Method from HAPPE v1 (NOT RECOMMENDED)', ...
                    'Choose Option', ...
                    'default', 'legacy', 'default');
           switch params.lineNoise.cl
                case 'default'
                        params.lineNoise.cl = 0;
                case 'legacy'
                        params.lineNoise.cl = 1;
            end 
        else
%             fprintf(['Enter lower cutoff, in Hz, for the notch filter:\n' ...
%                 'Example:' sprintf(' %i', params.lineNoise.freq-1) '\n']) ;
%             params.lineNoise.low = input('> ') ;
              lineNoise_low = inputdlg('Enter lower cutoff, in Hz, for the notch filter: ');
              params.lineNoise.low = str2num(lineNoise_low{1});
%             fprintf(['Enter higher cutoff, in Hz, for the notch filter:\n' ...
%                 'Example:' sprintf(' %i', params.lineNoise.freq+1) '\n']) ;
%             params.lineNoise.high = input('> ') ;
              lineNoise_high = inputdlg('Enter higher cutoff, in Hz, for the notch filter: ');
              params.lineNoise.high = str2num(lineNoise_high{1});
        end
    end
    
    %% RESAMPLE
    if (~preExist || strcmpi(paramChoice, 'resampling')) && ~reprocessing
        params.downsample = determ_downsample() ;
    end

    %% FILTER
    if params.paradigm.ERP.on && (~preExist || strcmpi(paramChoice, 'filter'))
        params.filt.on = 1 ;
%         params.filt.lowpass = input(['Enter low-pass filter, ' ...
%             'in Hz:\nCommon low-pass filter is 30 - 45 Hz.\n> ']) ;
        filt_lowpass = inputdlg('Enter low-pass filter, in Hz: Common low-pass filter is 30 - 45 Hz.') ;
        params.filt.lowpass = str2num(filt_lowpass{1});
%         params.filt.highpass = input(['Enter high-pass filter,' ...
%             ' in Hz:\nCommon high-pass filter is 0.1 - 0.3 Hz.\n> ']) ;
        filt_highpass = inputdlg('Enter high-pass filter, in Hz: Common low-pass filter is 0.1 - 0.3 Hz.') ;
        params.filt.highpass = str2num(filt_highpass{1});
        if params.loadInfo.chanlocs.inc
%             fprintf(['Choose a filter:\n  fir = Hamming windowed sinc FIR' ...
%                 'filter (EEGLAB''s standard filter)\n  butter = IIR ' ...
%                 'butterworth filter (ERPLAB''s standard filter)\n']) ;
%             params.filt.butter = choose2('fir', 'butter') ;
            params.filt.butter = questdlg('Choose a filter: fir = Hamming windowed sinc FIR filter (EEGLAB''s standard filter), butter = IIR butterworth filter (ERPLAB''s standard filter)', ...
                    'Choose Option', ...
                    'fir', 'butter', 'fir');
            switch params.filt.butter
                case 'fir'
                        params.filt.butter = 0;
                case 'butter'
                        params.filt.butter = 1;
            end 
        else; params.filt.butter = 0 ;
        end
    elseif ~params.paradigm.ERP.on && (~preExist || strcmpi(paramChoice, ...
            'filter') && ~reprocessing)
%         fprintf('Filter your data? [Y/N]\n') ;
%         params.filt.on = choose2('N', 'Y') ;
        params.filt.on = questdlg('Filter your data?', ...
                'Choose Option', ...
                'yes', 'no', 'yes');
        switch params.filt.on
            case 'no'
                    params.filt.on = 0;
            case 'yes'
                    params.filt.on = 1;
        end 
        if params.filt.on
            params.filt.lowpass = [];
            params.filt.highpass = [];
%             fprintf(['Filter type:\n  lowpass = Use a low-pass filter\n  ' ...
%                 'highpass = Use a high-pass filter\n  bandpass = Use ' ...
%                 'a band-pass filter\n']) ;
%             while true
%                 ui = input('> ', 's') ;
%                 if strcmpi(ui, 'lowpass')
%                     params.filt.method = [1,0] ; break ;
%                 elseif strcmpi(ui, 'highpass')
%                     params.filt.method = [0,1] ; break ;
%                 elseif strcmpi(ui, 'bandpass')
%                     params.filt.method = [1,1] ; break ;
%                 else; fprintf(['Invalid input: please enter "lowpass", ' ...
%                         '"highpass", or "bandpass" (without quotations).']) ;
%                 end
%             end
            params.filt.method = questdlg('Filter type:', ...
                'Choose Option', ...
                'lowpass', 'highpass', 'bandpass','bandpass');
            switch params.filt.method
                case 'lowpass'
                        params.filt.method = [1,0];
                case 'highpass'
                        params.filt.method = [0,1];
                case 'bandpass'
                        params.filt.method = [1,1];
            end 
            if params.filt.method(1)
%                 params.filt.lowpass = input(['Enter low-pass filter, in Hz:\n' ...
%                     'Common low-pass filter is 100 Hz.\n> ']) ;
                filt_lowpass = inputdlg('Enter low-pass filter, in Hz') ;
                params.filt.lowpass = str2num(filt_lowpass{1});
            end
            if params.filt.method(2)
%                 params.filt.highpass = input(['Enter high-pass filter, in Hz:\n' ...
%                     'Common high-pass filter is 1 Hz. NOTE: Use a >=1 Hz filter' ...
%                     ' if running muscIL.\n> ']) ;
                filt_highpass = inputdlg('Enter high-pass filter, in Hz') ;
                params.filt.highpass = str2num(filt_highpass{1});
            end
        else; params.filt.method = [0,0] ;
        end
        params.filt.butter = 0 ;
    end

    %% BAD CHANNEL DETECTION
    if (~preExist || strcmpi(paramChoice, 'bad channel detection')) && ~reprocessing
%         fprintf('Perform bad channel detection? [Y/N]\n') ;
%         params.badChans.rej = choose2('n','y') ;
        params.badChans.rej = questdlg('Perform bad channel detection?', ...
                'Choose Option', ...
                'yes', 'no', 'yes');
        switch params.badChans.rej
            case 'no'
                    params.badChans.rej = 0;
            case 'yes'
                    params.badChans.rej = 1;
        end 
        if params.badChans.rej
            if params.lowDensity; params.badChans.order = 0 ;
            else
%                 fprintf(['Run bad channel detection before or after wavelet' ...
%                     ' thresholding?\n  before = Run bad channel detection ' ...
%                     'first.\n  after = Run wavelet thresholding first (retains ' ...
%                     'more channels)\n']) ;
%                 params.badChans.order = choose2('before', 'after') ;
                params.badChans.order = questdlg('Run bad channel detection before or after wavelet thresholding?', ...
                'Choose Option', ...
                'before', 'after', 'after');
                switch params.badChans.order
                    case 'before'
                            params.badChans.order = 0;
                    case 'after'
                            params.badChans.order = 1;
                end 
            end
%             fprintf(['Run bad channel detection before or after wavelet' ...
%                 ' thresholding?\n  before = Run bad channel detection ' ...
%                 'first.\n  after = Run wavelet thresholding first (retains ' ...
%                 'more channels)\n']) ;
%             params.badChans.order = choose2('before', 'after') ;
        end
    end

    %% ECGONE
    if (~preExist || strcmpi(paramChoice, 'ECGone')) && ~reprocessing
%         fprintf(['Use ECGone (adapted from Isler et al., 2022) to reduce ' ...
%             'excess ECG artifact in your data? [Y/N]\n']) ;
%         params.ecgone.on = choose2('N', 'Y') ;
        params.ecgone.on = questdlg('Use ECGone (adapted from Isler et al., 2022) to reduce excess ECG artifact in your data?', ...
        'Choose Option', ...
        'yes', 'no', 'no');
        switch params.ecgone.on
            case 'no'
                    params.ecgone.on = 0;
            case 'yes'
                    params.ecgone.on = 1;
        end 
        if params.ecgone.on
%             fprintf(['Does your data contain a dedicated ECG channel? [Y/N]' ...
%                 '\n']) ;
%             params.ecgone.ECGchan.inc = choose2('N', 'Y') ;
            params.ecgone.ECGchan.inc = questdlg('Does your data contain a dedicated ECG channel?', ...
            'Choose Option', ...
            'yes', 'no', 'yes');
            switch params.ecgone.ECGchan.inc
                case 'no'
                        params.ecgone.ECGchan.inc = 0;
                case 'yes'
                        params.ecgone.ECGchan.inc = 1;
            end 
            if params.ecgone.ECGchan.inc
%                 fprintf('Dedicated ECG channel name:\n') ;
%                 params.ecgone.ECGchan.ID = input('> ', 's') ;
                ecgchanID = inputdlg('Dedicated ECG channel name:');
                params.ecgone.ECGchan.ID = ecgchanID{1};
                params.ecgone.peaky = 0 ;
            else
%                 fprintf(['Enter the channel names containing ECG artifact' ...
%                     ' to create a proxy:\nPress enter/return between each' ...
%                     ' entry.\nExample: E17\nWhen you have entered all ' ...
%                     'channels, input "done" (without quotations).\n']) ;
%                 params.ecgone.ECGchan.ID = unique(UI_cellArray(1, {}), 'stable') ;
                ecgchanID = inputdlg('Enter the channel names containing ECG artifact to create a proxy (separated by blank space):');
                params.ecgone.ECGchan.ID = unique(strsplit(ecgchanID{1}), 'stable');

%                 fprintf(['Threshold for determining if proxy contains ' ...
%                     'significant ECG artifact:\nExample: 10\n']) ;
%                 params.ecgone.peaky = input('> ') ;
                ecpeaky = inputdlg('Threshold for determining if proxy contains significant ECG artifact:');
                params.ecgone.peaky = str2num(ecpeaky{1});
            end
%             fprintf('Window creation length in MILLISECONDS for peak detection:\n') ;
%             params.ecgone.peakWinSize = input('> ')/1000 ;
            ecpeakWinSize = inputdlg('Window creation length in MILLISECONDS for peak detection:');
            params.ecgone.peakWinSize = str2num(ecpeakWinSize{1});
%             fprintf('Template window creation length in SECONDS:\n') ;
%             params.ecgone.procEpoch = input('> ') ;
            params.ecgone.procEpoch = 30 ;
        end
    end

    %% WAVELET METHODOLOGY
    if (~preExist || strcmpi(paramChoice, 'wavelet thresholding')) && ~reprocessing
%         fprintf(['Method of wavelet thresholding:\n  default = Default' ...
%             'method optimized in HAPPE v2.\n  legacy = Method from HAPPE' ...
%             ' v1 (NOT RECOMMENDED).\n']) ;
%         params.wavelet.legacy = choose2('default', 'legacy') ;
        params.wavelet.legacy = questdlg('Method of wavelet thresholding: default = Default method optimized in HAPPE v2, legacy = Method from HAPPE v1 (NOT RECOMMENDED)', ...
        'Choose Option', ...
        'default', 'legacy', 'default');
        switch params.wavelet.legacy
            case 'default'
                    params.wavelet.legacy = 0;
            case 'legacy'
                    params.wavelet.legacy = 1;
        end 
        if ~params.wavelet.legacy
            if params.paradigm.ERP.on
%                 fprintf(['Threshold rule for wavelet thresholding:\n  soft - ' ...
%                     'Use a soft threshold\n' ...
%                     '  hard - Use a hard threshold\n']) ;
%                 params.wavelet.softThresh = choose2('hard', 'soft') ;
                params.wavelet.softThresh = questdlg('Threshold rule for wavelet thresholding:', ...
                'Choose Option', ...
                'soft', 'hard', 'hard');
                switch params.wavelet.softThresh
                    case 'hard'
                            params.wavelet.softThresh = 0;
                    case 'soft'
                            params.wavelet.softThresh = 1;
                end 
            end
        end
    end
        
    %% MUSCIL
    if ~preExist || strcmpi(paramChoice, 'muscIL')
        if ~params.paradigm.ERP.on
%             fprintf(['Use ICLabel to reduce remaining muscle artifact ' ...
%                 'in your data? [Y/N]\nNOTE: This will drastically increase' ...
%                 ' processing time. Recommended for files with significant' ...
%                 ' muscle artifact.\n']) ;
%             params.muscIL = choose2('N', 'Y') ;
            params.muscIL = questdlg('Use ICLabel to reduce remaining muscle artifact in your data?', ...
            'Choose Option', ...
            'yes', 'no', 'no');
            switch params.muscIL
                case 'no'
                        params.muscIL = 0;
                case 'yes'
                        params.muscIL = 1;
            end 
        else; params.muscIL = 0;
        end
    end
    
    %% SEGMENTATION
    if ~preExist || strcmpi(paramChoice, 'segmentation')
%         fprintf('Segment data? [Y/N]\n') ;
%         params.segment.on = choose2('N', 'Y') ;
        params.segment.on = questdlg('Segment data?', ...
        'Choose Option', ...
        'yes', 'no', 'yes');
        switch params.segment.on
            case 'no'
                    params.segment.on = 0;
            case 'yes'
                    params.segment.on = 1;
        end 
        if params.segment.on
            if params.paradigm.task
%                 params.segment.start = input(['Segment start, in MILLISECONDS, ' ...
%                     'relative to stimulus onset:\nExample: -100\n> '])/1000 ;
%                 params.segment.end = input(['Segment end, in MILLISECONDS, ' ...
%                     'relative to stimulus onset:\n> '])/1000 ;
                segment_start = inputdlg('Segment start, in MILLISECONDS, relative to stimulus onset:') ;
                params.segment.start = str2num(segment_start{1})/1000;
                segment_end = inputdlg('Segment end, in MILLISECONDS, relative to stimulus onset:') ;
                params.segment.end = str2num(segment_end{1})/1000;
%                 params.segment.bounds = NaN(size(params.paradigm.onsetTags,2),2) ;
                % SET SEGMENT START AND END
%                 for i=1:size(params.segment.bounds,1)
%                     params.segment.bounds(1,1) = input(['Segment start, in MILLISECONDS, ' ...
%                         'relative to stimulus onset for "' params.paradigm.onsetTags{i} ...
%                         '":\nExample: -100\n> '])/1000 ;
%                     params.segment.bounds(1,2) = input(['Segment end, in MILLISECONDS, ' ...
%                         'relative to stimulus onset for "' params.paradigm.onsetTags{i} ...
%                         '":\n> '])/1000 ;
%                 end
                if params.paradigm.ERP.on
                    % DETERMINE TASK OFFSET
                    % *** For this, maybe make it possible to upload a list
                    % of offset delays?
%                     params.segment.offset = input(['Offset delay, in MILLISECONDS, ' ...
%                         'between stimulus initiation and presentation:\n' ...
%                         'NOTE: Please enter the total offset (combined system' ...
%                         ' and task-specific offsets).\n' ...
%                         '> ']) ;
                    segment_offset = inputdlg('Offset delay, in MILLISECONDS, between stimulus initiation and presentation:') ;
                    params.segment.offset = str2num(segment_offset{1});
                    % DETERMINE IF WANT BASELINE CORRECTION
%                     fprintf('Perform baseline correction (by subtraction)? [Y/N]\n') ;
%                     params.baseCorr.on = choose2('n', 'y') ;
                    params.baseCorr.on = questdlg('Perform baseline correction (by subtraction)?', ...
                    'Choose Option', ...
                    'yes', 'no', 'yes');
                    switch params.baseCorr.on
                        case 'no'
                                params.baseCorr.on = 0;
                        case 'yes'
                                params.baseCorr.on = 1;
                    end 

                    if params.baseCorr.on
                        % DETERMINE BASELINE START AND END
%                         params.baseCorr.start = input(['Enter, in MILLISECONDS,' ...
%                             ' where the baseline segment begins:\nExample: -100\n> ']) ;
                        baseCorr_start = inputdlg('Enter, in MILLISECONDS, where the baseline segment begins:') ;
                        params.baseCorr.start = str2num(baseCorr_start{1});
%                         params.baseCorr.end = input(['Enter, in MILLISECONDS,' ...
%                             ' where the baseline segment ends:\n' ...
%                             'NOTE: 0 indicates stimulus onset.\n> ']) ;
                        baseCorr_end = inputdlg('Enter, in MILLISECONDS, where the baseline segment ends:') ;
                        params.baseCorr.end = str2num(baseCorr_end{1});
                    end
                end
            % DETERMINE SEGMENT LENGTH
            elseif ~params.paradigm.task
%                 params.segment.length = input("Segment length, in SECONDS:\n> ") ;
                segment_length = inputdlg('Segment length, in SECONDS:') ;
                params.segment.length = str2num(segment_length{1});
    
            end
        end
    end
    
    %% INTERPOLATION
    if ~preExist || strcmpi(paramChoice, 'interpolation')
        if ~params.loadInfo.chanlocs.inc  % || ~params.badChans.rej
            params.segment.interp = 0 ;
        else
%             fprintf(['Interpolate the specific channels data determined ' ...
%                 'to be artifact/bad within each segment? [Y/N]\n']) ;
%             params.segment.interp = choose2('n', 'y') ;
            params.segment.interp = questdlg('Interpolate the specific channels data determined to be artifact/bad within each segment?', ...
            'Choose Option', ...
            'yes', 'no', 'yes');
            switch params.segment.interp
                case 'no'
                        params.segment.interp = 0;
                case 'yes'
                        params.segment.interp = 1;
            end 

        end
    end
    
    %% SEGMENT REJECTION
    if ~preExist || strcmpi(paramChoice, 'segment rejection')
%         fprintf('Perform segment rejection? [Y/N]\n') ;
%         params.segRej.on = choose2('n', 'y') ;
        params.segRej.on = questdlg('Perform segment rejection?', ...
        'Choose Option', ...
        'yes', 'no', 'yes');
        switch params.segRej.on
            case 'no'
                    params.segRej.on = 0;
            case 'yes'
                    params.segRej.on = 1;
        end 
        if params.segRej.on
%             fprintf(['Choose a method of segment rejection:\n  amplitude =' ...
%                 ' Amplitude criteria only\n  similarity = Segment ' ...
%                 'similarity only\n  both = Both amplitude criteria and ' ...
%                 'segment similarity\n']) ;
             params.segRej.method = questdlg('Choose a method of segment rejection: amplitude = Amplitude criteria only, similarity = Segment similarity only, both = Both amplitude criteria and segment similarity', ...
            'Choose Option', ...
            'amplitude', 'similarity', 'both', 'amplitude');
%             while true
%                 params.segRej.method = input('> ','s') ;
                if strcmpi(params.segRej.method, 'amplitude') || ...
                        strcmpi(params.segRej.method, 'both')
%                     params.segRej.minAmp = input(['Minimum signal amplitude' ...
%                         ' to use as the artifact threshold:\n> ']) ;
%                     params.segRej.maxAmp = input(['Maximum signal amplitude' ...
%                         'to use as the artifact threshold:\n> ']) ;
                    segRej_minAmp = inputdlg('Minimum signal amplitude to use as the artifact threshold:') ;
                    params.segRej.minAmp = str2num(segRej_minAmp{1});
                    segRej_maxAmp = inputdlg('Maximum signal amplitude to use as the artifact threshold:') ;
                    params.segRej.maxAmp = str2num(segRej_maxAmp{1});
%                     break ;
%                 elseif strcmpi(params.segRej.method, 'similarity'); break ;
%                 else
%                     fprintf(['Invalid input: please enter "amplitude", ' ...
%                         '"similarity", or "both" (without quotations)\n']) ;
                end
%             end
%             fprintf(['Use all channels or a region of interest for ' ...
%                 'segment rejection?\n  all = All channels\n  roi = Region' ...
%                 ' of interest\n']) ;
%             params.segRej.ROI.on = choose2('all', 'roi') ;
            params.segRej.ROI.on = questdlg('Use all channels or a region of interest for segment rejection? all = All channels, roi = Region of interest', ...
            'Choose Option', ...
            'all', 'roi', 'all');
            switch params.segRej.ROI.on
                case 'all'
                        params.segRej.ROI.on = 0;
                case 'roi'
                        params.segRej.ROI.on = 1;
            end 
            if params.segRej.ROI.on
%                 fprintf(['Choose an option for entering channels:\n  include = ' ...
%                     'Include ONLY the entered channel names\n  exclude = Include ' ...
%                     'every channel EXCEPT the entered channel names\n']) ;
%                 params.segRej.ROI.include = choose2('exclude', 'include') ;
                params.segRej.ROI.include = questdlg('Choose an option for entering channels: include = Include ONLY the entered channel names, exclude = Include every channel EXCEPT the entered channel names', ...
                'Choose Option', ...
                'include', 'exclude', 'include');
                switch params.segRej.ROI.include 
                    case 'exclude'
                            params.segRej.ROI.include = 0;
                    case 'include'
                            params.segRej.ROI.include = 1;
                end 
%                 fprintf(['Enter channels, including the preceding letter, one at ' ...
%                     'a time.\nPress enter/return between each entry.\nExamples: ' ...
%                     'E17\n          M1\nWhen you have entered all channels, input ' ...
%                     '"done" (without quotations).\n']) ;
%                 params.segRej.ROI.chans = unique(UI_cellArray(1, {}), 'stable') ;
                segRej_ROIchans = inputdlg('Enter channels, including the preceding letter, separated by a blank space') ;
%                 segRej_ROIchans = str2num(segRej_ROIchans{1});
                params.segRej.ROI.chans = unique(strsplit(segRej_ROIchans{1}), 'stable');

            end
            if params.paradigm.task && params.loadInfo.inputFormat == 1
%                 fprintf(['Use pre-selected "usable" trials to restrict ' ...
%                     'analysis? [Y/N]\n']) ;
%                 params.segRej.selTrials = choose2('n', 'y') ;
                params.segRej.selTrials = questdlg('Use pre-selected "usable" trials to restrict analysis?', ...
                'Choose Option', ...
                'yes', 'no', 'yes');
                switch params.segRej.selTrials 
                    case 'no'
                            params.segRej.selTrials = 0;
                    case 'yes'
                            params.segRej.selTrials = 1;
                end 
%                 if params.segRej.selTrials
%                     fprintf(['Enter the file, including the full path name ' ...
%                         'and file extension, indicating which trials should' ...
%                         ' be included in analyses.\n']) ;
%                     params.segRej.trialFile = input('> ', 's') ;
%                 end
            end
        end
    end
    
    %% RE-REFERENCING
    if ~preExist || strcmpi(paramChoice, 're-referencing')
        if ~params.loadInfo.chanlocs.inc
            params.reref.on = 0;
        else
%             fprintf('Re-reference data? [Y/N]\n') ;
%             params.reref.on = choose2('n', 'y') ;
            params.reref.on = questdlg('Re-reference data?', ...
            'Choose Option', ...
            'yes', 'no', 'yes');
            switch params.reref.on 
                case 'no'
                        params.reref.on = 0;
                case 'yes'
                        params.reref.on = 1;
            end 
            if params.reref.on
%                 fprintf(['Does your data contain a flatline or all zero ' ...
%                     'reference channel? [Y/N]\n']) ;
%                 params.reref.flat = choose2('n', 'y') ;
                params.reref.flat = questdlg('Does your data contain a flatline or all zero reference channel?', ...
                'Choose Option', ...
                'yes', 'no', 'yes');
                switch params.reref.flat
                    case 'no'
                            params.reref.flat = 0;
                    case 'yes'
                            params.reref.flat = 1;
                end 
                if params.reref.flat
%                     fprintf(['Enter reference channel ID:\nIf unknown, ' ...
%                         'press enter/return.\n']) ;
%                     params.reref.chan = input('> ', 's') ;
                     reref_chan = inputdlg('Enter reference channel ID:') ;
                     params.reref.chan = reref_chan{1};
                end
                
%                 fprintf(['Re-referencing type:\n  average = Average ' ...
%                     're-referencing\n  subset = Re-reference to a channel ' ...
%                     'or subset of channels\n  rest = Re-reference using' ...
%                     ' infinity with REST (Yao, 2001)\n']) ;
                params.reref.method = questdlg('Re-referencing type: average = Average  re-referencing, subset = Re-reference to a channel or subset of channels, rest = Re-reference using infinity with REST (Yao, 2001)', ...
                'Choose Option', ...
                'average', 'subset', 'rest', 'average');
%                 while true
%                     params.reref.method = input('> ', 's') ;
                    if strcmpi(params.reref.method, 'subset')
%                         fprintf(['Enter channel/subset of channels to ' ...
%                         're-reference to, one at a time.\nWhen you have ' ...
%                         'entered all channels, input "done" (without ' ...
%                         'quotations).\n']) ;
%                         params.reref.subset = UI_cellArray(1,{}) ;
                         reref_subset = inputdlg('Enter channel/subset of channels to re-reference to, separated by a blank space') ;
                         params.reref.subset = unique(strsplit(reref_subset{1}), 'stable');

%                         break ;
%                     elseif strcmpi(params.reref.method, 'average'); break;
%                     elseif strcmpi(params.reref.method, 'rest')
%                         fprintf(['Reference type:\n  average = Average ' ...
%                             'reference\n  fixed = Fixed reference\n  ' ...
%                             'mastoid = Ipsilateral and/or contralateral ' ...
%                             'mastoid\n']) ;
%                         while true
%                             params.reref.rest.ref = input('> ', 's') ;
%                             if strcmpi(params.reref.rest.ref, 'average') || ...
%                                 strcmpi(params.reref.rest.ref, 'fixed')
%                                 params.reref.rest.chanlist = {} ;
%                                 % Chanlist should be full list? Set as
%                                 % origchans
%                                 break ;
%                             elseif strcmpi(params.reref.rest.ref, 'mastoid')
%                                 while true
%                                     fprintf(['Enter list of left channels, ' ...
%                                         'one at a time.\nWhen ' ...
%                                         'finished entering channels, enter' ...
%                                         ' "done" (without quotations).\n']) ;
%                                     params.reref.rest.chanL = UI_cellArray(1, {});
%                         
%                                     fprintf(['Enter list of right channels, ' ...
%                                         'one at a time.\nWhen ' ...
%                                         'finished entering channels, enter' ...
%                                         ' "done" (without quotations).\n']) ;
%                                     params.reref.rest.chanR = UI_cellArray(1, {});
%                         
%                                     if isempty(intersect(params.reref.rest.chanL, ...
%                                             params.reref.rest.chanR)); break;
%                                     else
%                                         fprintf(['Invalid input: left and ' ...
%                                             'right channel lists cannot ' ...
%                                             'share any elements.\n']) ;
%                                     end
%                                 end
%                                 break ;
%                             else
%                                 fprintf(['Invalid input: please enter ' ...
%                                     '"average", "fixed", or "mastoid" ' ...
%                                     '(without quotations).\n']) ;
%                             end
%                         end
% 
%                         fprintf(['Load or calculate leadfield?\n  load = Load' ...
%                             ' existing leadfield\n  calculate = Calculate new' ...
%                             ' leadfield\n']) ;
%                         params.reref.rest.calc = choose2('load', 'calculate') ;
%                         if params.reref.rest.calc && strcmpi(params.reref.rest.ref, 'mastoid')
%                             while true
%                                 fprintf(['Enter XYZ coordinates for the' ...
%                                     ' left reference:\nEnter in this ' ...
%                                     'format: [x y z] (including brackets),' ...
%                                     ' then press enter/return.\n']) ;
%                                 params.reref.rest.coordLeft = input('> ') ;
%                                 if size(params.reref.rest.coordLeft, 2) == 3
%                                     break ;
%                                 else
%                                     fprintf(['Invalid input: please ' ...
%                                         'enter three numbers.\n']) ;
%                                 end
%                             end
%                             
%                             while true
%                                 fprintf(['Enter XYZ coordinates for the' ...
%                                     ' right reference:\nEnter in this ' ...
%                                     'format: [x y z] (including brackets),' ...
%                                     ' then press enter/return.\n']) ;
%                                 params.reref.rest.coordRight = input('> ') ;
%                                 if size(params.reref.rest.coordRight, 2) == 3
%                                     break ;
%                                 else
%                                     fprintf(['Invalid input: please ' ...
%                                         'enter three numbers.\n']) ;
%                                 end
%                             end
%                         elseif ~params.reref.rest.calc
%                             fprintf(['Enter the name of the leadfield file ' ...
%                                 'including path and file extension:\n']) ;
%                             while true
%                                 params.reref.rest.file = input('> ', 's') ;
%                                 if isfile(params.reref.rest.file); break;
%                                 else; fprintf(['Invalid input: please ' ...
%                                         'enter an existing file.\n']);
%                                 end
%                             end
%                         end
%                         break ;
%                     else
%                         fprintf(['Invalid input: please enter "average", ' ...
%                             '"subset", or "infinity" (without quotations)\n']) ;
                    end
            end
        end
    end
        break
    end
    
    %% SAVE FORMAT
    if ~preExist || strcmpi(paramChoice, 'save format')
        params.outputFormat = determ_saveFormat() ;
    end
    
    %% VISUALIZATIONS
    if ~preExist || strcmpi(paramChoice, 'visualizations')
%         fprintf('Run HAPPE with visualizations? [Y/N]\n') ;
%         params.vis.enabled = choose2('N', 'Y') ;
        params.vis.enabled = questdlg('Run HAPPE with visualizations?', ...
        'Choose Option', ...
        'yes', 'no', 'yes');
        switch params.vis.enabled
            case 'no'
                    params.segRej.vis.enabled = 0;
            case 'yes'
                    params.segRej.vis.enabled = 1;
        end 
        if params.vis.enabled
            % POWER SPECTRUM:
            % Min and Max
            vis_min = inputdlg("Minimum value for power spectrum figure:") ;
            params.vis.min = str2num(vis_min{1});
            vis_max = inputdlg("Maximum value for power spectrum figure:") ;
            params.vis.max = str2num(vis_max{1});
            % Frequencies for spatial topoplots
%             fprintf(['Enter the frequencies, one at a time, to generate ' ...
%                 'spatial topoplots for:\nWhen you have entered all ' ...
%                 'frequencies, input "done" (without quotations).\n']) ;
            vis_toPlot = inputdlg('Enter the frequencies, separated by a blank space, to generate spatial topoplots for:') ;
            params.vis.toPlot = unique(strsplit(vis_toPlot), 'stable');
%             indx = 1 ;
%             while true
%                 user_input = input('> ', 's') ;
%                 if strcmpi(user_input, 'done')
%                     if ~isfield(params.vis, 'toPlot')
%                          params.vis.toPlot = [];
%                     end
%                     params.vis.toPlot = unique(params.vis.toPlot, 'stable') ;
%                     break ;
%                 else
%                     params.vis.toPlot(indx) = str2double(user_input) ;
%                     indx = indx + 1 ;
%                 end
%             end
            
            if params.paradigm.ERP.on
                % DETERMINE TIME RANGE FOR THE TIMESERIES FIGURE        
%                 params.vis.min = input('Start time, in MILLISECONDS, for the ERP timeseries figure:\n> ') ;
%                 params.vis.max = input(['End time, in MILLISECONDS, for the ERP timeseries figure:\n' ...
%                     'NOTE: This should end 1+ millisecond(s) before your segmentation parameter ends (e.g. 299 for 300).\n' ...
%                     '> ']) ;
                vis_min = inputdlg('Start time, in MILLISECONDS, for the ERP timeseries figure:') ;
                params.vis.min = str2num(vis_min{1});
                vis_max = inputdlg('End time, in MILLISECONDS, for the ERP timeseries figure: NOTE: This should end 1+ millisecond(s) before your segmentation parameter ends (e.g. 299 for 300).') ;
                params.vis.max = str2num(vis_max{1});
                % Frequencies for spatial topoplots
%                 fprintf(['Enter the latencies, one at a time, to ' ...
%                     'generate spatial topoplots for:\nWhen you have ' ...
%                     'entered all latencies, input "done" (without ' ...
%                     'quotations).\n']) ;
%                 indx = 1 ;
%                 while true
%                     user_input = input('> ', 's') ;
%                     if strcmpi(user_input, 'done')
%                         params.vis.toPlot = unique(params.vis.toPlot, 'stable') ;
%                         break ;
%                     else
%                         params.vis.toPlot(indx) = str2double(user_input) ;
%                         indx = indx + 1 ;
%                     end
%                 end
                vis_toPlot = inputdlg('Enter the latencies, separated by a blank space, to generate spatial topoplots for:') ;
                params.vis.toPlot = unique(strsplit(vis_toPlot), 'stable');
            end
        end
    end
    
   %% DONE
   if ~preExist || strcmpi(paramChoice, 'done')
%        fprintf('Please check your parameters before continuing.\n') ;
%        listParams(params) ;
        gui = makeParamGUI();      
        listParamsGUI(gui, params)
%        fprintf('Are the above parameters correct? [Y/N]\n') ;
        params_correct = questdlg('Are the parameters correct?', ...
        'Choose Option', ...
        'yes', 'no', 'yes');
        switch params_correct
            case 'yes'
                    return ;
            case 'no'
                    if ~preExist
                        changedParams = 1 ;
                        preExist = 1 ;
                    end
        end 
%        if choose2('n','y'); break ;
%        elseif ~preExist
%            changedParams = 1 ;
%            preExist = 1 ;
%        end
   end   
end