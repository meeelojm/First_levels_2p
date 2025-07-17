function fps = read_fps(folder,defaultFPS)
    ini = dir(fullfile(folder,'*.ini'));
    fps = defaultFPS;
    if ~isempty(ini)
        txt = fileread(fullfile(folder,ini(1).name));
        tok = regexp(txt,'volume\.rate\.\(in\.Hz\)\s*=\s*([0-9.]+)','tokens','once');
        if ~isempty(tok), fps = str2double(tok{1}); end
    end
end