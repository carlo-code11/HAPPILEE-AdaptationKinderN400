function guiPrint(gui, txt)
%helper function to print to a gui
    old = gui.out.Value;
    if ischar(old)
        old = {old};
    end
    gui.out.Value = [old; {txt}];
    drawnow;
end
