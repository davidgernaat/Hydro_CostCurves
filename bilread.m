function [X,Y ] = bilread(bilname, latlim,lonlim, rows,cols)
% Read data from BIL (band interleaved) raster file.
% Optionally, only a small window, identified by either LATLIM/LONLIM or ROWS/COLS may be read.

% Copyright P.W.Bogaart 2015,2016
% bilname = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\data\EUR\bil\3s\bilcon\n50e025_con_bil\n50e025_con.bil';
% latlim
% lonlim
% rows
% cols

bilname;
isfile = exist(bilname, 'file');

if isfile==0
    X = 0;
    Y = 0;
    disp('No 3s file')
    bilname
    return;
end

% read all raster data, or just a block selection?
if nargin==5
    assert(isempty(latlim));
    assert(isempty(lonlim));
    partial = 2;
elseif nargin==3
    partial = 1;
else
    partial = 0;
end;

% Read header
hdr = bilheader(bilname);

nrows = hdr.nrows;
ncols = hdr.ncols;
if partial>0
%     fprintf('Reading BIL file %s (%dx%d) [partially]\n', bilname, nrows,ncols);
else
%     fprintf('Reading BIL file %s (%dx%d) [full]\n', bilname, nrows,ncols);
end


% determine wordsize (i.e. number of bytes per BIL element)
if     hdr.nbits ==32, wordsize=4; filefmt='int32'; iofmt='*int32';
elseif hdr.nbits ==16, wordsize=2; filefmt='int16'; iofmt='*int16';
elseif hdr.nbits == 8, wordsize=1; filefmt='uint8'; iofmt='*uint8';
else
    error(sprintf('Invalid nbits=%d'));
end;

if partial>0
    if partial==1 % Use lat/lon limits to specifity range
        
        latlim;
        lonlim;
        
        % lat,lon limits should be a simple range
        assert(numel(latlim)==2);
        assert(numel(lonlim)==2);
        
        % enforce range to be increasing: [min, max]
        if latlim(2)<latlim(1),
            latlim = [latlim(2), latlim(1)];
        end;
        if lonlim(2)<lonlim(1),
            lonlim = [lonlim(2), lonlim(1)];
        end;
        
        % Coordinates of tile boundaries
        west  = hdr.ulxmap - 0.5*hdr.xdim;
        east  = west + hdr.ncols * hdr.xdim;
        north = hdr.ulymap + 0.5*hdr.ydim;
        south = north - hdr.nrows * hdr.ydim;
        
        % see if the requested window are aligned with cell boundaries or cell
        % centers
        r1 = (hdr.ydim + hdr.ulymap - latlim(2)) / hdr.ydim;
        r2 = (hdr.ydim + hdr.ulymap - latlim(1)) / hdr.ydim;
        c1 = (lonlim(1) + hdr.xdim - hdr.ulxmap) ./ hdr.xdim;
        c2 = (lonlim(2) + hdr.xdim - hdr.ulxmap) ./ hdr.xdim;
        if approx(r1-fix(r1), 0.5)
            aligned = 'B';
            firstrow = round(r1+0.5);
            lastrow  = round(r2-0.5);
            firstcol = round(c1+0.5);
            lastcol  = round(c2-0.5);
        elseif approx(r1-fix(r1), 0.0)
            aligned = 'C';
            firstrow = round(r1);
            lastrow  = round(r2);
            firstcol = round(c1);
            lastcol  = round(c2);
        else
            error('Invalid lims');
        end
        
    elseif partial==2 % use rows/cols to specify range
        
        firstrow = rows(1);
        lastrow  = rows(end);
        firstcol = cols(1);
        lastcol  = cols(end);
    end;
    
    % Check if requested rows/cols are valid, given DEM size
    assert(firstrow>=1);
    assert(firstcol>=1);
    assert(lastrow<=nrows);
    assert(lastcol<=ncols);
    
    % How much rows/cols to read or skip?
    nrowread  = lastrow - firstrow +1;
    nrowskip = firstrow-1;
    
    ncolread = lastcol - firstcol+1;
    ncolskip1 = firstcol-1;
    ncolskip2 = ncols-lastcol;
    
    % load DEM, region of interest (ROI) only
    dem = zeros(nrowread, ncolread,filefmt);
    
    fid = fopen(bilname,'r');
    fseek(fid, nrowskip*ncols*wordsize, 'bof'); % skip over rows above ROI
    for r = 1:nrowread
        fseek(fid, ncolskip1*wordsize, 0); % skip columns left of ROI
        row = fread(fid, ncolread, iofmt); % read columns within ROI
        dem(r,:) = row;
        fseek(fid, ncolskip2*wordsize, 0); % skip columns right of ROI
    end;
    fclose(fid);
    
else % Real ALL data
    
    fprintf('Loading: %s (%dx%d) [%s]...', bilname, nrows,ncols, iofmt);
    fid = fopen(bilname, 'rb');
    dem = fread(fid, [ncols nrows],iofmt)'; % read rows as columns, then transpose...
    fclose(fid);
    fprintf('\n');
    
end

X = dem;
if nargout>1
    Y = hdr;
end
end

function y = snap(x, target, info)
if approx(x, target)
    fprintf('Snapping to %s\n', info);
    y = target;
else
    y = x;
end
end

function good = approx(x, target)
good = abs(x-target) < 0.001*target;
end