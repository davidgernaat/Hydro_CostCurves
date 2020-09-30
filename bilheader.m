function hdr = bilheader(bilname)

[path, name, ~] = fileparts(bilname);

% Read header
hdrname = fullfile(path, [name, '.hdr']);
assert(exist(hdrname,'file')==2);

fid = fopen(hdrname, 'rt');
C = textscan(fid,'%s %s');
fclose(fid);

hdr.fname = hdrname;

n = length(C{1});
for i=1:n
    key = lower(C{1}{i});
    val = C{2}{i};
    if     strcmp(key,'nrows'),  hdr.nrows  = str2double(val);
    elseif strcmp(key,'ncols'),  hdr.ncols  = str2double(val);
    elseif strcmp(key,'ulxmap'), hdr.ulxmap = str2double(val);
    elseif strcmp(key,'ulymap'), hdr.ulymap = str2double(val);
    elseif strcmp(key,'xdim'),   hdr.xdim   = str2double(val);
    elseif strcmp(key,'ydim'),   hdr.ydim   = str2double(val);
    elseif strcmp(key,'nbits'),  hdr.nbits  = str2double(val);
    elseif strcmp(key,'nodata'), hdr.nodata = str2double(val);
    else
        % fprintf('Ignoring %s = %s\n', key, val);
    end
end

%% Make xdim and ydim more precise
hdr.xdim = snap(hdr.xdim,  3/3600, '3 arc-second');
hdr.xdim = snap(hdr.xdim, 15/3600, '15 arc-second');
hdr.xdim = snap(hdr.xdim, 30/3600, '30 arc-second');
hdr.xdim = snap(hdr.xdim, 1/60, '1 arc-minute');
hdr.xdim = snap(hdr.xdim, 5/60, '5 arc-minute');

hdr.ydim = snap(hdr.ydim,  3/3600, '3 arc-second');
hdr.ydim = snap(hdr.ydim, 15/3600, '15 arc-second');
hdr.ydim = snap(hdr.ydim, 30/3600, '30 arc-second');
hdr.ydim = snap(hdr.ydim, 1/60, '1 arc-minute');
hdr.ydim = snap(hdr.ydim, 5/60, '5 arc-minute');

end

function y = snap(x, target, info)
if approx(x, target)
    %fprintf('Snapping to %s\n', info);
    y = target;
else
    y = x;
end
end

function good = approx(x, target)
good = abs(x-target) < 0.001*target;
end