function gui = makeParamGUI()
%helper function to create a GUI
    gui.fig = uifigure('Name','Parameter Settings','Position',[300 100 500 710]);
    gui.out = uitextarea(gui.fig,...
        'Position',[10 10 480 700],...
        'Editable','off');
end


    
