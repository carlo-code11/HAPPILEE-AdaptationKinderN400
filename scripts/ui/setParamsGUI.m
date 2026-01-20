function params = setParamsGUI(params, preExist, reprocessing, changedParams, happeDir)
% setParamsGUI  Programmatic GUI for editing HAPPE parameters
% Usage: params = setParamsGUI(params, preExist, reprocessing, changedParams, happeDir)
%
% Opens a uifigure GUI allowing the user to set/modify parameters using
% buttons, dropdowns, checkboxes instead of command-line input. Returns the
% updated params when the user clicks Confirm, or the original params if
% the user clicks Cancel or closes the window.

arguments
    params struct = struct()
    preExist (1,1) {mustBeNumericOrLogical} = 0
    reprocessing (1,1) {mustBeNumericOrLogical} = 0
    changedParams (1,1) {mustBeNumericOrLogical} = 0
    happeDir char = ''
end

% Ensure some fields exist with sensible defaults
params = ensureDefaults(params);

% Create UIFigure
fig = uifigure('Name','HAPPE - Set Parameters','Position',[200 100 1100 700]);
% Use a grid layout: left controls, right summary
grid = uigridlayout(fig,[1,2]);
grid.ColumnWidth = {'0.58x','0.42x'};
grid.RowHeight = {'1x'};

% Left panel: controls (scrollable)
leftPanel = uipanel(grid,'Title','Parameters','Padding',[8 8 8 8]);
leftPanel.Layout.Row = 1;
leftPanel.Layout.Column = 1;

scroll = uiflowcontainer(leftPanel,'FlowDirection','topdown','Overflow','on');
scroll.FlowDirection = 'topdown';
scroll.Width = leftPanel.Position(3)-16;
scroll.Height = leftPanel.Position(4)-16;

% Right panel: summary + advanced editor + buttons
rightPanel = uipanel(grid,'Title','Summary & Actions','Padding',[8 8 8 8]);
rightPanel.Layout.Row = 1;
rightPanel.Layout.Column = 2;

% Build sections on left
% --- Density
lbl = uilabel(scroll,'Text','Density (channels)','FontWeight','bold');
densityGroup = uibuttongroup(scroll,'Title','Low density (<=30) ?','Padding',[8 8 8 8]);
densityGroup.FlowBreak = 'off';
densityGroup.Layout.Width = scroll.Width - 20;
rbLow = uiradiobutton(densityGroup,'Text','Low (<=30)','Position',[10 10 120 22]);
rbHigh = uiradiobutton(densityGroup,'Text','High (>30)','Position',[140 10 120 22]);
if isfield(params,'lowDensity') && params.lowDensity
    densityGroup.SelectedObject = rbLow;
else
    densityGroup.SelectedObject = rbHigh;
end
densityGroup.SelectionChangedFcn = @onDensityChange;

% --- Paradigm (rest/task) + ERP
lbl = uilabel(scroll,'Text','Paradigm','FontWeight','bold');
paradigmGrid = uigridlayout(scroll,[1,3]);
paradigmGrid.ColumnWidth = {'1x','100px','100px'};
paradigmGrid.Padding = [0 0 0 0];
paradigmGrid.Layout.Width = scroll.Width - 20;

ddParadigm = uidropdown(paradigmGrid,'Items',{'rest','task'});
if isfield(params,'paradigm') && isfield(params.paradigm,'task') && params.paradigm.task
    ddParadigm.Value = 'task';
else
    ddParadigm.Value = 'rest';
end
ddParadigm.ValueChangedFcn = @onParadigmChange;

lblERP = uilabel(paradigmGrid,'Text','ERP:');
tglERP = uiswitch(paradigmGrid,'on','Items',{'Off','On'});
if isfield(params,'paradigm') && isfield(params.paradigm,'ERP') && isfield(params.paradigm.ERP,'on') && params.paradigm.ERP.on
    tglERP.Value = 'On';
else
    tglERP.Value = 'Off';
end
tglERP.ValueChangedFcn = @onERPChange;

% --- Onset tags (if task)
lbl = uilabel(scroll,'Text','Task Onset Tags','FontWeight','bold');
hBox = uigridlayout(scroll,[1,2]);
hBox.ColumnWidth = {'1x','120px'};
listOnset = uilistbox(hBox,'Items',{params.paradigm.onsetTags{:}},'Multiselect','off');
listOnset.Layout.Row = 1; listOnset.Layout.Column = 1;
panelButtons = uipanel(hBox,'Title','','Padding',[4 4 4 4]);
panelButtons.Layout.Row = 1; panelButtons.Layout.Column = 2;
btnAddOnset = uibutton(panelButtons,'Text','Add','ButtonPushedFcn',@onAddOnset);
btnRemOnset = uibutton(panelButtons,'Text','Remove','ButtonPushedFcn',@onRemoveOnset,'Position',[10 50 100 22]);
edtOnset = uieditfield(panelButtons,'text','Position',[10 85 100 22],'Placeholder','tag name');

% --- Channels of Interest
lbl = uilabel(scroll,'Text','Channels of Interest','FontWeight','bold');
chanGrid = uigridlayout(scroll,[1,3]);
chanGrid.ColumnWidth = {'1x','1x','1x'};
btnChAll = uibutton(chanGrid,'Text','All','ButtonPushedFcn',@(~,~) setChannels('all'));
btnChInclude = uibutton(chanGrid,'Text','Include...','ButtonPushedFcn',@(~,~) setChannels('include'));
btnChExclude = uibutton(chanGrid,'Text','Exclude...','ButtonPushedFcn',@(~,~) setChannels('exclude'));

% --- Line noise
lbl = uilabel(scroll,'Text','Line Noise / Harmonics','FontWeight','bold');
lnGrid = uigridlayout(scroll,[2,3]);
lnGrid.ColumnWidth = {'1x','1x','1x'};
lblFreq = uilabel(lnGrid,'Text','Line frequency (Hz):'); edtFreq = uieditfield(lnGrid,'numeric','Value',getFieldOr(params,'lineNoise.freq',60));
lblFreq.Layout.Row = 1; lblFreq.Layout.Column = 1; edtFreq.Layout.Row = 1; edtFreq.Layout.Column = 2;
lblHarm = uilabel(lnGrid,'Text','Harmonics (comma list):'); edtHarm = uieditfield(lnGrid,'text','Value',strjoin(arrayfun(@num2str,getFieldOr(params,'lineNoise.harms.freqs',[]),'UniformOutput',false),','));
lblMethod = uilabel(lnGrid,'Text','Method:'); ddLN = uidropdown(lnGrid,'Items',{'notch','cleanline'},'Value',getFieldOr(params,'lineNoise.cl',false,'mapBoolToMethod'));
ddLN.Layout.Row = 2; ddLN.Layout.Column = 3;

% --- Resampling
lbl = uilabel(scroll,'Text','Resample','FontWeight','bold');
resampleGrid = uigridlayout(scroll,[1,3]);
ddResample = uidropdown(resampleGrid,'Items',{'None','To custom (Hz)'},'Value', getFieldOr(params,'downsample',0,'mapDownsample'));
edtResample = uieditfield(resampleGrid,'numeric','Value',getFieldOr(params,'downsample',0));
ddResample.ValueChangedFcn = @onResampleChange;

% --- Filtering
lbl = uilabel(scroll,'Text','Filtering','FontWeight','bold');
filtGrid = uigridlayout(scroll,[2,3]);
filtGrid.ColumnWidth = {'1x','1x','1x'};
lblFiltOn = uilabel(filtGrid,'Text','Filter On?'); ddFiltOn = uidropdown(filtGrid,'Items',{'No','Yes'},'Value',ternary(getFieldOr(params,'filt.on',0)==1,'Yes','No'));
lblFiltOn.Layout.Row = 1; lblFiltOn.Layout.Column = 1;
lblLow = uilabel(filtGrid,'Text','Lowpass (Hz):'); edtLow = uieditfield(filtGrid,'numeric','Value',getFieldOr(params,'filt.lowpass',[]));
edtLow.Layout.Row=2; edtLow.Layout.Column=1;
lblHigh = uilabel(filtGrid,'Text','Highpass (Hz):'); edtHigh = uieditfield(filtGrid,'numeric','Value',getFieldOr(params,'filt.highpass',[]));
edtHigh.Layout.Row=2; edtHigh.Layout.Column=2;
ddFiltOn.ValueChangedFcn = @onFiltChange;

% --- Bad channel detection
lbl = uilabel(scroll,'Text','Bad Channel Detection','FontWeight','bold');
badGrid = uigridlayout(scroll,[1,2]);
btnBadOn = uibutton(badGrid,'Text','Toggle On/Off','ButtonPushedFcn',@onBadToggle);
lblBadState = uilabel(badGrid,'Text',ternary(getFieldOr(params,'badChans.rej',0),'On','Off'));

% --- ECGone
lbl = uilabel(scroll,'Text','ECGone','FontWeight','bold');
ecgGrid = uigridlayout(scroll,[2,2]);
btnECG = uibutton(ecgGrid,'Text','Enable/Disable','ButtonPushedFcn',@onECGToggle);
lblECGstate = uilabel(ecgGrid,'Text',ternary(getFieldOr(params,'ecgone.on',0),'On','Off'));
btnECG.Layout.Row=1; btnECG.Layout.Column=1; lblECGstate.Layout.Row=1; lblECGstate.Layout.Column=2;

% --- Wavelet
lbl = uilabel(scroll,'Text','Wavelet Thresholding','FontWeight','bold');
waveGrid = uigridlayout(scroll,[1,3]);
ddWave = uidropdown(waveGrid,'Items',{'default','legacy'},'Value', ternary(getFieldOr(params,'wavelet.legacy',1)==1,'legacy','default'));
ddWave.ValueChangedFcn = @onWaveChange;

% --- muscIL
lbl = uilabel(scroll,'Text','muscIL (muscle ICA)','FontWeight','bold');
muscBtn = uibutton(scroll,'Text',ternary(getFieldOr(params,'muscIL',0),'Disable','Enable'),'ButtonPushedFcn',@onMuscToggle);

% --- Segmentation
lbl = uilabel(scroll,'Text','Segmentation','FontWeight','bold');
segGrid = uigridlayout(scroll,[2,3]);
ddSegOn = uidropdown(segGrid,'Items',{'No','Yes'},'Value',ternary(getFieldOr(params,'segment.on',0)==1,'Yes','No'));
edtSegStart = uieditfield(segGrid,'numeric','Value',getFieldOr(params,'segment.start',-0.1));
edtSegEnd = uieditfield(segGrid,'numeric','Value',getFieldOr(params,'segment.end',0.6));
ddSegOn.ValueChangedFcn = @onSegToggle;

% --- Interpolation
lbl = uilabel(scroll,'Text','Interpolation','FontWeight','bold');
ddInterp = uidropdown(scroll,'Items',{'No','Yes'},'Value', ternary(getFieldOr(params,'segment.interp',0)==1,'Yes','No'));
ddInterp.ValueChangedFcn = @(~,~) setField(params,'segment.interp',strcmp(ddInterp.Value,'Yes'));

% --- Segment Rejection (method)
lbl = uilabel(scroll,'Text','Segment Rejection','FontWeight','bold');
rejGrid = uigridlayout(scroll,[1,3]);
ddRejOn = uidropdown(rejGrid,'Items',{'No','Yes'},'Value',ternary(getFieldOr(params,'segRej.on',0)==1,'Yes','No'));
ddRejMethod = uidropdown(rejGrid,'Items',{'off','amplitude','similarity','both'},'Value', ternary(getFieldOr(params,'segRej.method','off'),'off'));
ddRejOn.ValueChangedFcn = @onRejToggle;

% --- Re-Referencing
lbl = uilabel(scroll,'Text','Re-Referencing','FontWeight','bold');
ddRerefOn = uidropdown(scroll,'Items',{'No','Yes'},'Value',ternary(getFieldOr(params,'reref.on',0)==1,'Yes','No'));
ddRerefMethod = uidropdown(scroll,'Items',{'average','subset','rest'},'Value',getFieldOr(params,'reref.method','average'));
ddRerefOn.ValueChangedFcn = @onRerefToggle;

% --- Visualizations
lbl = uilabel(scroll,'Text','Visualizations','FontWeight','bold');
ddVis = uidropdown(scroll,'Items',{'No','Yes'},'Value',ternary(getFieldOr(params,'vis.enabled',0)==1,'Yes','No'));
edtVisMin = uieditfield(scroll,'numeric','Value',getFieldOr(params,'vis.min',[]));
edtVisMax = uieditfield(scroll,'numeric','Value',getFieldOr(params,'vis.max',[]));

% --- Save format
lbl = uilabel(scroll,'Text','Save Format','FontWeight','bold');
ddSave = uidropdown(scroll,'Items',{'txt','.mat','.set'},'Value',ternary(getFieldOr(params,'outputFormat',1)==1,'.txt',ternary(getFieldOr(params,'outputFormat',2)==2,'.mat','.set')));

% Add some spacing filler
spacer = uilabel(scroll,'Text','');

% Right side: summary & advanced editor & action buttons
txtArea = uitextarea(rightPanel,'Position',[10 120 rightPanel.Position(3)-20 rightPanel.Position(4)-140],'Editable','off');
txtArea.FontName = 'Courier New';
btnRefresh = uibutton(rightPanel,'Text','Refresh Summary','Position',[10 80 120 28],'ButtonPushedFcn',@refreshSummary);
btnAdvanced = uibutton(rightPanel,'Text','Advanced Editor','Position',[140 80 120 28],'ButtonPushedFcn',@openAdvancedEditor);

btnApply = uibutton(rightPanel,'Text','Apply','Position',[10 30 80 36],'ButtonPushedFcn',@onApply);
btnConfirm = uibutton(rightPanel,'Text','Confirm','Position',[100 30 80 36],'ButtonPushedFcn',@onConfirm);
btnCancel = uibutton(rightPanel,'Text','Cancel','Position',[190 30 80 36],'ButtonPushedFcn',@onCancel);

% Initialize display
refreshSummary();

% Block until Confirm/Cancel
uiwait(fig);

% After resume, get the final params from fig.UserData
if isvalid(fig) && isfield(fig,'UserData') && ~isempty(fig.UserData)
    params = fig.UserData;
else
    % If canceled or closed without setting UserData, return original params
    % (params already holds original unless modified)
end
if isvalid(fig)
    close(fig);
end

%% -------------------- Nested helper callbacks and functions --------------------

    function refreshSummary(~,~)
        % update summary text area to reflect current params
        try
            txtArea.Value = strsplit(renderParamsText(params), '\n');
        catch ME
            txtArea.Value = {'Error rendering summary:', ME.message};
        end
    end

    function onDensityChange(~,~)
        params.lowDensity = strcmp(densityGroup.SelectedObject.Text,'Low (<=30)');
        refreshSummary();
    end

    function onParadigmChange(~,~)
        params.paradigm.task = strcmp(ddParadigm.Value,'task');
        % if switching to rest, ensure ERP off
        if ~params.paradigm.task
            params.paradigm.ERP.on = 0;
            tglERP.Value = 'Off';
        end
        refreshSummary();
    end

    function onERPChange(~,~)
        params.paradigm.ERP.on = strcmp(tglERP.Value,'On');
        refreshSummary();
    end

    function onAddOnset(~,~)
        newTag = strtrim(edtOnset.Value);
        if isempty(newTag)
            uialert(fig,'Please type an onset tag before clicking Add.','Missing input');
            return;
        end
        if ~isfield(params.paradigm,'onsetTags') || isempty(params.paradigm.onsetTags)
            params.paradigm.onsetTags = {newTag};
        else
            params.paradigm.onsetTags{end+1} = newTag;
        end
        listOnset.Items = params.paradigm.onsetTags;
        edtOnset.Value = '';
        refreshSummary();
    end

    function onRemoveOnset(~,~)
        sel = listOnset.Value;
        if isempty(sel), return; end
        idx = find(strcmp(params.paradigm.onsetTags,sel),1);
        if ~isempty(idx)
            params.paradigm.onsetTags(idx) = [];
            listOnset.Items = params.paradigm.onsetTags;
        end
        refreshSummary();
    end

    function setChannels(mode)
        switch lower(mode)
            case 'all'
                params.chans.subset = 'all';
                params.chans.IDs = {};
            case 'include'
                % prompt a list dialog to add channels
                c = inputdlg('Enter channel IDs (comma-separated):','Include channels',1,{strjoin(getFieldOr(params,'chans.IDs',{}),',')});
                if ~isempty(c)
                    ids = strtrim(split(c{1},','));
                    params.chans.subset = 'coi_include';
                    params.chans.IDs = ids(:)';
                end
            case 'exclude'
                c = inputdlg('Enter channel IDs to exclude (comma-separated):','Exclude channels',1,{strjoin(getFieldOr(params,'chans.IDs',{}),',')});
                if ~isempty(c)
                    ids = strtrim(split(c{1},','));
                    params.chans.subset = 'coi_exclude';
                    params.chans.IDs = ids(:)';
                end
        end
        refreshSummary();
    end

    function onResampleChange(~,~)
        if strcmp(ddResample.Value,'None')
            params.downsample = 0;
            edtResample.Value = 0;
            edtResample.Editable = 'off';
        else
            edtResample.Editable = 'on';
            params.downsample = edtResample.Value;
        end
        refreshSummary();
    end

    function onFiltChange(~,~)
        params.filt.on = strcmp(ddFiltOn.Value,'Yes');
        params.filt.lowpass = edtLow.Value;
        params.filt.highpass = edtHigh.Value;
        refreshSummary();
    end

    function onBadToggle(~,~)
        current = getFieldOr(params,'badChans.rej',0);
        params.badChans.rej = ~current;
        lblBadState.Text = ternary(params.badChans.rej,'On','Off');
        refreshSummary();
    end

    function onECGToggle(~,~)
        params.ecgone.on = ~getFieldOr(params,'ecgone.on',0);
        lblECGstate.Text = ternary(params.ecgone.on,'On','Off');
        refreshSummary();
    end

    function onWaveChange(~,~)
        params.wavelet.legacy = strcmp(ddWave.Value,'legacy');
        refreshSummary();
    end

    function onMuscToggle(src,~)
        params.muscIL = ~getFieldOr(params,'muscIL',0);
        src.Text = ternary(params.muscIL,'Disable','Enable');
        refreshSummary();
    end

    function onSegToggle(~,~)
        params.segment.on = strcmp(ddSegOn.Value,'Yes');
        params.segment.start = edtSegStart.Value;
        params.segment.end = edtSegEnd.Value;
        refreshSummary();
    end

    function onRejToggle(~,~)
        params.segRej.on = strcmp(ddRejOn.Value,'Yes');
        params.segRej.method = ddRejMethod.Value;
        refreshSummary();
    end

    function onRerefToggle(~,~)
        params.reref.on = strcmp(ddRerefOn.Value,'Yes');
        params.reref.method = ddRerefMethod.Value;
        refreshSummary();
    end

    function onApply(~,~)
        % Apply current GUI values into params (synchronize)
        params.lineNoise.freq = edtFreq.Value;
        if isempty(edtHarm.Value)
            params.lineNoise.harms.freqs = [];
            params.lineNoise.harms.on = 0;
        else
            freqs = strtrim(split(edtHarm.Value,','));
            params.lineNoise.harms.freqs = cellfun(@str2double, freqs);
            params.lineNoise.harms.on = ~isempty(params.lineNoise.harms.freqs);
        end
        params.lineNoise.cl = strcmp(ddLN.Value,'cleanline');
        % resample
        if strcmp(ddResample.Value,'None')
            params.downsample = 0;
        else
            params.downsample = edtResample.Value;
        end
        % filter
        params.filt.on = strcmp(ddFiltOn.Value,'Yes');
        params.filt.lowpass = edtLow.Value;
        params.filt.highpass = edtHigh.Value;
        % segmentation
        params.segment.on = strcmp(ddSegOn.Value,'Yes');
        params.segment.start = edtSegStart.Value;
        params.segment.end = edtSegEnd.Value;
        % interpolation
        params.segment.interp = strcmp(ddInterp.Value,'Yes');
        % re-referencing
        params.reref.on = strcmp(ddRerefOn.Value,'Yes');
        params.reref.method = ddRerefMethod.Value;
        % vis
        params.vis.enabled = strcmp(ddVis.Value,'Yes');
        params.vis.min = edtVisMin.Value;
        params.vis.max = edtVisMax.Value;
        % save
        switch ddSave.Value
            case '.txt', params.outputFormat = 1;
            case '.mat', params.outputFormat = 2;
            case '.set', params.outputFormat = 3;
        end
        refreshSummary();
        uialert(fig,'Changes applied to parameters (not yet confirmed).','Applied');
    end

    function onConfirm(~,~)
        % Apply final fields then set fig.UserData and resume
        onApply();
        fig.UserData = params;
        uiresume(fig);
    end

    function onCancel(~,~)
        % Close without saving (return original params)
        fig.UserData = [];
        uiresume(fig);
    end

    function openAdvancedEditor(~,~)
        % Open a modal dialog to edit the entire params struct as JSON.
        d = uifigure('Name','Advanced Params Editor','Position',[300 150 600 500]);
        ta = uitextarea(d,'Position',[10 50 580 430],'FontName','Courier New');
        try
            ta.Value = jsonencode(params);
        catch
            ta.Value = sprintf('Error encoding params to JSON\n');
        end
        btnOK = uibutton(d,'Text','OK','Position',[10 10 80 30],'ButtonPushedFcn',@(~,~) advOK());
        btnCancelAdv = uibutton(d,'Text','Cancel','Position',[100 10 80 30],'ButtonPushedFcn',@(~,~) close(d));
        function advOK()
            try
                newparams = jsondecode(ta.Value);
                % convert struct-like (if needed)
                params = newparams;
                close(d);
                refreshSummary();
            catch ME
                uialert(d, ['Error parsing JSON: ' ME.message],'JSON Error');
            end
        end
    end

%% -------------------- Utility functions --------------------

    function out = getFieldOr(s, fld, default, varargin)
        % getFieldOr(s,'a.b.c',default) safe nested get
        if nargin < 3, default = []; end
        parts = split(fld,'.');
        cur = s;
        ok = true;
        for k = 1:numel(parts)
            p = parts{k};
            if isstruct(cur) && isfield(cur,p)
                cur = cur.(p);
            else
                ok = false; break;
            end
        end
        if ok, out = cur; else out = default; end
        % optional mapping modes
        if nargin>3 && strcmp(varargin{1},'mapBoolToMethod')
            if out, out = 'cleanline'; else out = 'notch'; end
        elseif nargin>3 && strcmp(varargin{1},'mapDownsample')
            if isempty(out) || out==0, out = 'None'; else out = 'To custom (Hz)'; end
        end
    end

    function s = ternary(cond, a, b)
        if cond, s = a; else s = b; end
    end

    function params = ensureDefaults(params)
        % populate minimal fields so GUI won't error
        if ~isfield(params,'lowDensity'), params.lowDensity = 0; end
        if ~isfield(params,'paradigm'), params.paradigm = struct(); end
        if ~isfield(params.paradigm,'task'), params.paradigm.task = 0; end
        if ~isfield(params.paradigm,'ERP'), params.paradigm.ERP.on = 0; end
        if ~isfield(params.paradigm,'onsetTags'), params.paradigm.onsetTags = {}; end
        if ~isfield(params,'QCfreqs'), params.QCfreqs = [2,5,8,12,20,30,45,70]; end
        if ~isfield(params,'loadInfo'), params.loadInfo = struct('inputFormat',1,'chanlocs',struct('inc',0)); end
        if ~isfield(params,'chans'), params.chans = struct('subset','all','IDs',{{}}); end
        if ~isfield(params,'lineNoise'), params.lineNoise = struct('freq',60,'harms',struct('on',0,'freqs',[]),'cl',0); end
        if ~isfield(params,'downsample'), params.downsample = 0; end
        if ~isfield(params,'filt'), params.filt = struct('on',0,'lowpass',[],'highpass',[],'butter',0); end
        if ~isfield(params,'badChans'), params.badChans = struct('rej',0,'order',0); end
        if ~isfield(params,'ecgone'), params.ecgone = struct('on',0,'ECGchan',struct('inc',0,'ID',{}),'peaky',[],'procEpoch',30,'peakWinSize',1); end
        if ~isfield(params,'wavelet'), params.wavelet = struct('legacy',1,'softThresh',0); end
        if ~isfield(params,'muscIL'), params.muscIL = 0; end
        if ~isfield(params,'segment'), params.segment = struct('on',0,'start',-0.1,'end',0.6,'interp',0,'length',1); end
        if ~isfield(params,'segRej'), params.segRej = struct('on',0,'method','off','minAmp',[],'maxAmp',[],'ROI',struct('on',0,'include',1,'chans',{{}}),'selTrials',0); end
        if ~isfield(params,'reref'), params.reref = struct('on',0,'method','average','flat',0,'chan',''); end
        if ~isfield(params,'vis'), params.vis = struct('enabled',0,'min',[],'max',[],'toPlot',[]); end
        if ~isfield(params,'outputFormat'), params.outputFormat = 1; end
    end

    function txt = renderParamsText(p)
        % A readable text representation of params (condensed version of listParams)
        % This function is intentionally simpler than the original long print.
        lines = {};
        lines{end+1} = '---------------------------------------------';
        lines{end+1} = 'PARAMETER SETTINGS:';
        % Density
        lines{end+1} = sprintf('Density: %s', ternary(getFieldOr(p,'lowDensity',0),'Low (<= 30 channels)','High (> 30 channels)'));
        % Paradigm
        if getFieldOr(p,'paradigm.task',0)
            lines{end+1} = 'Paradigm: Task';
            lines{end+1} = [' - Onset Tags: ' strjoin(getFieldOr(p,'paradigm.onsetTags',{}),', ')];
            lines{end+1} = sprintf(' - ERP Analysis: %s', ternary(getFieldOr(p,'paradigm.ERP.on',0),'Yes','No'));
        else
            lines{end+1} = 'Paradigm: Resting State';
        end
        % Data format
        fmt = getFieldOr(p,'loadInfo.inputFormat',1);
        fmtName = {' .mat (NetStation/MATLAB)',' .raw (Netstation)',' .set (EEGLAB)',' .cdt (Neuroscan)',' .mff (EGI)',' .edf',' .bdf'};
        idx = min(max(1,fmt),numel(fmtName));
        lines{end+1} = ['Data File Format:' fmtName{idx}];
        % Acquisition layout
        L = getFieldOr(p,'loadInfo.layout',[5, NaN]);
        if isnumeric(L) && numel(L)>=2
            switch L(1)
                case 1, lines{end+1} = sprintf('Acquisition: %i channel EGI Geodesic Sensor Net', L(2));
                case 2, lines{end+1} = sprintf('Acquisition: %i channel EGI HydroCel', L(2));
                case 3, lines{end+1} = sprintf('Acquisition: %i channel Neuroscan Quik-Cap', L(2));
                otherwise, lines{end+1} = 'Acquisition: Unspecified';
            end
        end
        % Channels
        cs = getFieldOr(p,'chans.subset','all');
        if strcmpi(cs,'all')
            lines{end+1} = 'Channels: All';
        elseif strcmpi(cs,'coi_include')
            lines{end+1} = ['Channels (include): ' strjoin(getFieldOr(p,'chans.IDs',{}),', ')];
        else
            lines{end+1} = ['Channels (exclude): ' strjoin(getFieldOr(p,'chans.IDs',{}),', ')];
        end
        % Line noise
        if isfield(p,'lineNoise')
            lines{end+1} = sprintf('Line Noise Frequency: %i Hz', getFieldOr(p,'lineNoise.freq',60));
            if getFieldOr(p,'lineNoise.harms.on',0)
                lines{end+1} = [' - Additional Frequencies: ' strjoin(arrayfun(@num2str,getFieldOr(p,'lineNoise.harms.freqs',[]),'UniformOutput',false),', ')];
            end
            lines{end+1} = sprintf('Line Noise Method: %s', ternary(getFieldOr(p,'lineNoise.cl',0),'CleanLine','Notch Filter'));
        end
        % Resample
        ds = getFieldOr(p,'downsample',0);
        lines{end+1} = sprintf('Resample: %s', ternary(ds==0,'Off',sprintf('To %i Hz',ds)));
        % Filter
        if getFieldOr(p,'filt.on',0)
            lines{end+1} = sprintf('Filter: On (Lowpass %s Hz, Highpass %s Hz)', num2str(getFieldOr(p,'filt.lowpass',NaN)), num2str(getFieldOr(p,'filt.highpass',NaN)));
        else
            lines{end+1} = 'Filter: Off';
        end
        % Bad channels
        lines{end+1} = sprintf('Bad Channel Detection: %s', ternary(getFieldOr(p,'badChans.rej',0),'On','Off'));
        % ECGone
        lines{end+1} = sprintf('ECGone: %s', ternary(getFieldOr(p,'ecgone.on',0),'On','Off'));
        % Wavelet
        lines{end+1} = sprintf('Wavelet: %s', ternary(getFieldOr(p,'wavelet.legacy',1),'Legacy','Default'));
        % MuscIL
        lines{end+1} = sprintf('MuscIL: %s', ternary(getFieldOr(p,'muscIL',0),'On','Off'));
        % Segmentation
        lines{end+1} = sprintf('Segmentation: %s', ternary(getFieldOr(p,'segment.on',0),'On','Off'));
        if getFieldOr(p,'segment.on',0)
            lines{end+1} = sprintf(' - Start: %g s  End: %g s', getFieldOr(p,'segment.start',NaN), getFieldOr(p,'segment.end',NaN));
            if getFieldOr(p,'paradigm.ERP.on',0)
                lines{end+1} = sprintf(' - Baseline correction: %s', ternary(getFieldOr(p,'baseCorr.on',0),'On','Off'));
            end
        end
        % Seg Rejection
        lines{end+1} = sprintf('Segment Rejection: %s', ternary(getFieldOr(p,'segRej.on',0),'On','Off'));
        % Re-referencing
        lines{end+1} = sprintf('Re-referencing: %s', ternary(getFieldOr(p,'reref.on',0),'On','Off'));
        % Visualizations
        lines{end+1} = sprintf('Visualizations: %s', ternary(getFieldOr(p,'vis.enabled',0),'On','Off'));
        % Save format
        fmt = getFieldOr(p,'outputFormat',1);
        sfmt = '.txt'; if fmt==2, sfmt='.mat'; elseif fmt==3, sfmt='.set'; end
        lines{end+1} = sprintf('Save Format: %s', sfmt);
        lines{end+1} = '---------------------------------------------';
        txt = strjoin(lines,'\n');
    end

end
