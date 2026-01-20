function listParamsGUI(gui, params)

guiPrint(gui, '---------------------------------------------');
guiPrint(gui, 'PARAMETER SETTINGS:');

%% Density
if params.lowDensity
    guiPrint(gui, 'Density: Low (<= 30 channels)');
else
    guiPrint(gui, 'Density: High (> 30 channels)');
end

%% Rest vs Task
if params.paradigm.task
    guiPrint(gui, 'Paradigm: Task');
    switch class(params.paradigm.onsetTags)
        case 'cell'
             guiPrint(gui, [' - Task Onset Tags: ' strjoin(params.paradigm.onsetTags, ', ')]);
        case 'char'
             guiPrint(gui, [' - Task Onset Tags: ' params.paradigm.onsetTags]);
    end
    guiPrint(gui, ' - Conditions:');
    if params.paradigm.conds.on
        for i = 1:size(params.paradigm.conds.groups,1)
            label = params.paradigm.conds.groups{i,1};
            items = params.paradigm.conds.groups(i,2:end);
            items = items(~cellfun('isempty',items));
            guiPrint(gui, ['   ' label ': ' strjoin(items, ', ')]);
        end
    else
        guiPrint(gui, 'NA');
    end
    
    if params.paradigm.ERP.on
        guiPrint(gui,' - ERP Analysis: Yes');
    else
        guiPrint(gui,' - ERP Analysis: No');
    end

else
    guiPrint(gui, 'Paradigm: Resting State');
end

%% Data File Format
guiPrint(gui, 'Data File Format:');
switch params.loadInfo.inputFormat
    case 1
        guiPrint(gui,'.mat (Net Station or MATLAB)');
        if params.loadInfo.NSformat
            guiPrint(gui,[' - EEG Variable Names: ' strjoin(params.loadInfo.NSvarNames, ', ')]);
        else
            if params.loadInfo.chanlocs.inc
                guiPrint(gui,[' - Channel Locations File: ' params.loadInfo.chanlocs.file]);
            else
                guiPrint(gui,' - Channel Locations Provided: No');
            end
            
            if params.loadInfo.srate.same
                guiPrint(gui,' - Same Sampling Rate: Yes');
            else
                guiPrint(gui,[' - Sampling Rates: ' params.loadInfo.srate.file]);
            end
        end
        
        if params.paradigm.task
            guiPrint(gui,[' - Task Event Info: ' params.loadInfo.eventLoc]);
        end
        
    case 2
        guiPrint(gui,'.raw (Netstation Simple Binary)');
    case 3
        guiPrint(gui,'.set (EEGLab)');
    case 4
        guiPrint(gui,'.cdt (Neuroscan)');
    case 5
        guiPrint(gui,'.mff (EGI)');
        guiPrint(gui,[' - Type Fields: ' strjoin(params.loadInfo.typeFields, ', ')]);
    case 6
        guiPrint(gui,'.edf');
    case 7
        guiPrint(gui,'.bdf → .set (Mentalab)');
end

%% Acquisition Layout
L = params.loadInfo.layout;
guiPrint(gui,'Acquisition Layout');
switch L(1)
    case 1
        guiPrint(gui, sprintf('%i-channel EGI Geodesic Sensor Net', L(2)));
    case 2
        guiPrint(gui, sprintf('%i-channel EGI HydroCel Net', L(2)));
    case 3
        guiPrint(gui, sprintf('%i-channel Neuroscan Quik-Cap', L(2)));
    case 4
        guiPrint(gui, sprintf('%i-channel Mentalab Explore', L(2)));
    otherwise
        guiPrint(gui,'Unspecified');
end

%% Channels
switch lower(params.chans.subset)
    case 'all'
        guiPrint(gui,'Channels: All');
    case 'coi_include'
        guiPrint(gui,['Channels: ' strjoin(params.chans.IDs, ', ')]);
    case 'coi_exclude'
        guiPrint(gui,['Channels: All except ' strjoin(params.chans.IDs, ', ')]);
end

%% Line Noise
if isfield(params,'lineNoise')
    guiPrint(gui, ['Line Noise Frequency: ' num2str(params.lineNoise.freq) ' Hz']);
    
    if params.lineNoise.harms.on
        freqs = arrayfun(@num2str, params.lineNoise.harms.freqs, 'UniformOutput', false);
        guiPrint(gui, [' - Additional Frequencies: ' strjoin(freqs, ', ')]);
    end
    
    if params.lineNoise.cl
        if params.lineNoise.legacy
            guiPrint(gui,'Line Noise Method: CleanLine (Legacy)');
        else
            guiPrint(gui,'Line Noise Method: CleanLine (Default)');
        end
    else
        guiPrint(gui,'Line Noise Method: Notch Filter');
        guiPrint(gui, [' - Low Cutoff: ' num2str(params.lineNoise.low)]);
        guiPrint(gui, [' - High Cutoff: ' num2str(params.lineNoise.high)]);
    end
end

%% Resampling
if isfield(params,'downsample')
    if params.downsample == 0
        guiPrint(gui,'Resample: Off');
    else
        guiPrint(gui,['Resample: to ' num2str(params.downsample) ' Hz']);
    end
end

%% Filtering
if params.filt.on
    guiPrint(gui,'Filter: On');
    guiPrint(gui,[' - Lowpass: ' num2str(params.filt.lowpass)]);
    guiPrint(gui,[' - Highpass: ' num2str(params.filt.highpass)]);
    if params.filt.butter
        guiPrint(gui,' - Type: ERPLAB Butterworth');
    else
        guiPrint(gui,' - Type: EEGLAB FIR');
    end
else
    guiPrint(gui,'Filter: Off');
end

%% Bad Channel Detection
if isfield(params,'badChans')
    if params.badChans.rej
        guiPrint(gui,'Bad Channel Detection: On');
        if params.badChans.order
            guiPrint(gui,' - Order: After Wavelet');
        else
            guiPrint(gui,' - Order: Before Wavelet');
        end
    else
        guiPrint(gui,'Bad Channel Detection: Off');
    end
end

%% ECGone
if isfield(params,'ecgone')
    if params.ecgone.on
        guiPrint(gui,'ECGone: On');
        if params.ecgone.ECGchan.inc
            guiPrint(gui,[' - ECG Channel: ' params.ecgone.ECGchan.ID]);
        else
            guiPrint(gui,[' - ECG Proxy Channels: ' strjoin(params.ecgone.ECGchan.ID, ', ')]);
            guiPrint(gui,[' - Artifact Threshold: ' num2str(params.ecgone.peaky)]);
        end
        guiPrint(gui,[' - Peak Window: ' num2str(params.ecgone.peakWinSize) ' s']);
        guiPrint(gui,[' - Template Length: ' num2str(params.ecgone.procEpoch) ' s']);
    else
        guiPrint(gui,'ECGone: Off');
    end
end

%% Wavelet
if isfield(params,'wavelet')
    if params.wavelet.legacy
        guiPrint(gui,'Wavelet Thresholding: Legacy');
    else
        guiPrint(gui,'Wavelet Thresholding: Default');
        if params.paradigm.ERP.on
            if params.wavelet.softThresh
                guiPrint(gui,' - Threshold Rule: Soft');
            else
                guiPrint(gui,' - Threshold Rule: Hard');
            end
        end
    end
end

%% MuscIL
if params.muscIL
    guiPrint(gui,'MuscIL: On');
else
    guiPrint(gui,'MuscIL: Off');
end

%% Segmentation
if params.segment.on
    guiPrint(gui,'Segmentation: On');
    
    if params.paradigm.task
        guiPrint(gui,[' - Start: ' num2str(params.segment.start) ' s']);
        guiPrint(gui,[' - End: ' num2str(params.segment.end) ' s']);
        
        if params.paradigm.ERP.on
            guiPrint(gui,[' - Offset: ' num2str(params.segment.offset) ' ms']);
            if params.baseCorr.on
                guiPrint(gui,' - Baseline Correction: On');
                guiPrint(gui,['   Start: ' num2str(params.baseCorr.start) ' ms']);
                guiPrint(gui,['   End: ' num2str(params.baseCorr.end) ' ms']);
            else
                guiPrint(gui,' - Baseline Correction: Off');
            end
        end
        
    else
        guiPrint(gui,[' - Segment Length: ' num2str(params.segment.length) ' s']);
    end
    
else
    guiPrint(gui,'Segmentation: Off');
end

%% Interpolation
if params.segment.interp
    guiPrint(gui,'Interpolation: On');
else
    guiPrint(gui,'Interpolation: Off');
end

%% Segment Rejection
if params.segRej.on
    guiPrint(gui,'Segment Rejection: On');
    guiPrint(gui,[' - Method: ' params.segRej.method]);
    
    if any(strcmpi(params.segRej.method, {'both','amplitude'}))
        guiPrint(gui,['   Min Amp: ' num2str(params.segRej.minAmp)]);
        guiPrint(gui,['   Max Amp: ' num2str(params.segRej.maxAmp)]);
        
        if params.segRej.ROI.on
            guiPrint(gui,'   ROI Mode');
            switch class(params.segRej.ROI.chans) 
                case 'cell'
                    guiPrint(gui,['   Channels: ' strjoin(params.segRej.ROI.chans, ', ')]);
                case 'char'
                    guiPrint(gui,['   Channels: ' params.segRej.ROI.chans]);
            end
        else
            guiPrint(gui,'   Using All Channels');
        end
        
        if params.paradigm.task && params.loadInfo.inputFormat == 1
            if params.segRej.selTrials
                guiPrint(gui,'   Selected Trials Only: Yes');
            else
                guiPrint(gui,'   Selected Trials Only: No');
            end
        end
    end
else
    guiPrint(gui,'Segment Rejection: Off');
end

%% Re-Referencing
if params.reref.on
    guiPrint(gui,'Re-Referencing: On');
    guiPrint(gui,[' - Method: ' params.reref.method]);
    
    switch lower(params.reref.method)
        case 'subset'
            guiPrint(gui,[' - Subset: ' strjoin(params.reref.subset, ', ')]);
            
        case 'rest'
            guiPrint(gui,[' - Reference Type: ' params.reref.rest.ref]);
            
            if strcmpi(params.reref.rest.ref,'mastoid')
                guiPrint(gui,['   Left: ' strjoin(params.reref.rest.chanL, ', ')]);
                guiPrint(gui,['   Right: ' strjoin(params.reref.rest.chanR, ', ')]);
            end
            
            if params.reref.rest.calc
                guiPrint(gui,' - Leadfield: Calculate');
            else
                guiPrint(gui,[' - Leadfield: Load (' params.reref.rest.file ')']);
            end
    end
    
else
    guiPrint(gui,'Re-Referencing: Off');
end

%% Visualizations
if params.vis.enabled
    guiPrint(gui,'Visualizations: On');
    guiPrint(gui,[' - Start: ' num2str(params.vis.min)]);
    guiPrint(gui,[' - End: ' num2str(params.vis.max)]);
    
    times = arrayfun(@num2str, params.vis.toPlot, 'UniformOutput', false);
    guiPrint(gui,[' - Times/Frequencies: ' strjoin(times, ', ')]);
else
    guiPrint(gui,'Visualizations: Off');
end

%% Save Format
switch params.outputFormat
    case 1
        guiPrint(gui,'Save Format: .txt');
    case 2
        guiPrint(gui,'Save Format: .mat');
    case 3
        guiPrint(gui,'Save Format: .set');
end

guiPrint(gui, '---------------------------------------------');

end
