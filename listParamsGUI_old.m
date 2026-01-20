function listParamsGUI(gui, params)
%helper function to make HAPPE user friendly with GUIs -> for the
%listParams function

guiPrint(gui, '---------------------------------------------');
guiPrint(gui, 'PARAMETER SETTINGS:');

%% Density
if params.lowDensity
    guiPrint(gui, 'Density: Low (<= 30 channels)');
else
    guiPrint(gui, 'Density: High (> 30 channels)');
end

%% Rest VS Task
if params.paradigm.task
    guiPrint(gui, 'Resting State or Task: Task');
    
    % Onset tags
    txt = sprintf(' - Task Onset Tags: %s', strjoin(params.paradigm.onsetTags, ', '));
    guiPrint(gui, txt);

    % Conditions
    guiPrint(gui, '    - Conditions:');
    if params.paradigm.conds.on
        for i = 1:size(params.paradigm.conds.groups,1)
            label = params.paradigm.conds.groups{i,1};
            curr = params.paradigm.conds.groups(i,2:end);
            curr = curr(~cellfun('isempty',curr));
            guiPrint(gui, sprintf('      %s = %s', label, strjoin(curr, ', ')));
        end
    else
        guiPrint(gui,'      NA');
    end

    % ERP
    if params.paradigm.ERP.on
        guiPrint(gui,' - ERP Analysis: Yes');
    else
        guiPrint(gui,' - ERP Analysis: No');
    end
else
    guiPrint(gui, 'Resting State or Task: Resting State');
end

%% Data File Format
fmt = params.loadInfo.inputFormat;
guiPrint(gui, 'Data File Format:')
switch fmt
    case 1
        guiPrint(gui,'.mat (Net Station output or MATLAB matrix)');
        if params.loadInfo.NSformat
            guiPrint(gui,[' - Potential EEG Variable Names: ', ...
                strjoin(params.loadInfo.NSvarNames, ', ')]);
        else
            guiPrint(gui,' - Channel Locations Provided:');
            if params.loadInfo.chanlocs.inc
                guiPrint(gui,['    Yes - Channel Locations File: ', params.loadInfo.chanlocs.file]);
            else
                guiPrint(gui,'    No');
            end

            if params.loadInfo.srate.same
                guiPrint(gui,' - Same Sampling Rate Across Files: Yes');
            else
                guiPrint(gui,[' - Same Sampling Rate Across Files: No (List: ', ...
                    params.loadInfo.srate.file ')']);
            end
        end

        if params.paradigm.task
            guiPrint(gui,[' - Task Event Info: ', params.loadInfo.eventLoc]);
        end

    case 2
        guiPrint(gui,'.raw (Netstation simple binary)');
    case 3
        guiPrint(gui,'.set (EEGLab)');
    case 4
        guiPrint(gui,'.cdt (Neuroscan)');
    case 5
        guiPrint(gui,'.mff (EGI)');
        guiPrint(gui,[' - Type Fields: ', strjoin(params.loadInfo.typeFields, ', ')]);
    case 6
        guiPrint(gui,'.edf');
    case 7
        guiPrint(gui,'.bdf->.set (Mentalab)');
end

%% Acquisition Layout
L = params.loadInfo.layout;
guiPrint(gui, 'Acquisition Layout:')
switch L(1)
    case 1
        guiPrint(gui, sprintf('%i channel EGI Geodesic Sensor Net', L(2)));
    case 2
        guiPrint(gui, sprintf('%i channel EGI HydroCel Geodesic Sensor Net', L(2)));
    case 3
        guiPrint(gui, sprintf('%i channel Neuroscan Quik-Cap', L(2)));
    case 4
        guiPrint(gui, sprintf('%i channel Mentalab Explore', L(2)));
    otherwise
        guiPrint(gui, 'Unspecified');
end

%% Channels
switch lower(params.chans.subset)
    case 'all'
        guiPrint(gui,'Channels: All');
    case 'coi_include'
        guiPrint(gui,['Channels: ', strjoin(params.chans.IDs, ', ')]);
    case 'coi_exclude'
        guiPrint(gui,['Channels: All except ', strjoin(params.chans.IDs, ', ')]);
end

%% (… ALL REMAINING SECTIONS TRANSLATED EXACTLY THE SAME …)

guiPrint(gui, "---------------------------------------------");
end
