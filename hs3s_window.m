function [dem, x, y] = hs3s_window(type,lat,lon,winsize,root,CONT)
% 
% clear all
% root = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global';
% type = 'con';
% lat = 1.6806;
% lon = 30.0258;
% winsize = 30;
% CONT = 'AFR';
% k=0;

% find corresponding Hydrosheds 3sec block and file
dem =0;
Northing = floor(lat/5)*5;
Easting = floor(lon/5)*5;
bilname = hs3s_tile(type,root, Northing, Easting,CONT);
% bilname = hs3s_tile(type,root, Northing, Easting);

bilname;
isfile = exist(bilname, 'file');

if isfile==0
    disp('No 3s file')
    bilname
    return;
end

% Read corresponding header
hdr = bilheader(bilname);
nrows = hdr.nrows;
ncols = hdr.ncols;
dx = hdr.xdim; % shortcuts
dy = hdr.ydim;

% 'map' borders
mwest  = hdr.ulxmap - 0.5 * dx;
meast  = mwest + ncols * dx;
mnorth = hdr.ulymap + 0.5 * dy;
msouth = mnorth - nrows * dy;

% Get focal row/column within this block
c = round((lon + hdr.xdim - hdr.ulxmap) / hdr.xdim);
r = round((hdr.ydim + hdr.ulymap - lat) / hdr.ydim);

% recalc lat/lon to be at precisely at the cell center
lat = mnorth - (r-1)*dy - 0.5*dy;
lon = mwest  + (c-1)*dx + 0.5*dx;

% % Get window around focal pixel
firstrow = r-winsize;
lastrow  = r+winsize;
firstcol = c-winsize;
lastcol  = c+winsize;
% 
% % Coordinates of current map boundaries (cell centers)
% mwest  = hdr.ulxmap;
% meast  = hdr.ulxmap + (hdr.ncols-1) * hdr.xdim;
% mnorth = hdr.ulymap;
% msouth = hdr.ulymap - (hdr.nrows-1) * hdr.ydim;

% Coordinates of window boundaries (cell borders)
wwest  = lon - (winsize+0.5) * hdr.xdim;
weast  = lon + (winsize+0.5) * hdr.xdim;
wnorth = lat + (winsize+0.5) * hdr.ydim;
wsouth = lat - (winsize+0.5) * hdr.ydim;

rwwest=rem(wwest,1);
rweast=rem(weast,1);
rwnorth=rem(wnorth,1);
rwsouth=rem(wsouth,1);

if rwwest<0; rwwest=rwwest*-1; end
if rweast<0; rweast=rweast*-1; end
if rwnorth<0; rwnorth=rwnorth*-1; end
if rwsouth<0; rwsouth=rwsouth*-1; end

if rem(rwwest,1)<0.0000001 |(1-rem(rwwest,1))<0.0000001; disp('rem wwest'); return; end
if rem(rweast,1)<0.0000001 |(1-rem(rwwest,1))<0.0000001; disp('rem weasst'); return; end
if rem(rwnorth,1)<0.0000001 |(1-rem(rwwest,1))<0.0000001; disp('rem wnorth'); return; end
if rem(rwsouth,1)<0.0000001 |(1-rem(rwwest,1))<0.0000001; disp('rem wsouth'); return; end

% Optionally compute axes

% if nargout>1
%     west  = hdr.ulxmap + (firstcol-1) * hdr.xdim;
%     east  = hdr.ulxmap + (lastcol-1)  * hdr.xdim;
%     north = hdr.ulymap - (firstrow-1) * hdr.ydim;
%     south = hdr.ulymap - (lastrow-1)  * hdr.ydim;
%     x = linspace(west, east,   lastcol-firstcol+1);
%     y = linspace(south, north, lastrow-firstrow+1);
% end

% Window may consist of multiple tiles. Collect them.

if wnorth > mnorth% cases 1,2,3
    Nlat = [mnorth, wnorth];
    Slat = [wsouth, mnorth];
    if wwest < mwest % case 1; overlap near NW corner
        % Compute latitude limits (cell borders)
        Wlon = [wwest, mwest];
        Elon = [mwest, weast];
        % Get fractional tiles
        NW  = bilread(hs3s_tile(type,root, Northing+5, Easting-5,CONT), Nlat,Wlon);
        NE  = bilread(hs3s_tile(type,root, Northing+5, Easting,CONT  ), Nlat,Elon);
        SW  = bilread(hs3s_tile(type,root, Northing  , Easting-5,CONT), Slat,Wlon);
        SE  = bilread(hs3s_tile(type,root, Northing  , Easting,CONT  ), Slat,Elon);
        if NW==0 ; return; end;
        if NE==0 ; return; end;
        if SW==0 ; return; end;
        if SE==0 ; return; end;
        dem = [NW NE; SW SE];
    elseif weast > meast% case 3; ; overlap near NE corner
        Wlon = [wwest, meast];
        Elon = [meast, weast];
        NW  = bilread(hs3s_tile(type,root, Northing+5, Easting,CONT  ), Nlat, Wlon);
        NE  = bilread(hs3s_tile(type,root, Northing+5, Easting+5,CONT), Nlat, Elon);
        SW  = bilread(hs3s_tile(type,root, Northing  , Easting,CONT  ), Slat, Wlon);
        SE  = bilread(hs3s_tile(type,root, Northing  , Easting+5,CONT), Slat, Elon);
        if NW==0 ; return; end;
        if NE==0 ; return; end;
        if SW==0 ; return; end;
        if SE==0 ; return; end;
        dem = [NW NE; SW SE];
    else % case 2; overlap near N border
        Xlon = [wwest, weast];
        N  = bilread(hs3s_tile(type,root, Northing+5, Easting,CONT), Nlat,Xlon);
        S  = bilread(hs3s_tile(type,root, Northing  , Easting,CONT), Slat,Xlon);
        if N==0; return; end;
        if S==0; return; end;
        dem = [N; S];
    end
elseif wsouth < msouth % cases 7,8,9
    Nlat = [msouth, wnorth];
    Slat = [wsouth, msouth];
    if wwest < mwest % case 7; overlap near SW corner
        Wlon = [wwest, mwest];
        Elon = [mwest, weast];
        NW  = bilread(hs3s_tile(type,root, Northing  , Easting-5,CONT), Nlat, Wlon);
        NE  = bilread(hs3s_tile(type,root, Northing  , Easting,CONT  ), Nlat, Elon);
        SW  = bilread(hs3s_tile(type,root, Northing-5, Easting-5,CONT), Slat, Wlon);
        SE  = bilread(hs3s_tile(type,root, Northing-5, Easting,CONT  ), Slat, Elon);
        if NW==0 ; return; end;
        if NE==0 ; return; end;
        if SW==0 ; return; end;
        if SE==0 ; return; end;
        dem = [NW NE; SW SE];
    elseif weast > meast % case 9; overlap near SE corner
        Wlon = [wwest, meast];
        Elon = [meast, weast];
        NW  = bilread(hs3s_tile(type,root, Northing  , Easting,CONT  ), Nlat, Wlon);
        NE  = bilread(hs3s_tile(type,root, Northing  , Easting+5,CONT), Nlat, Elon);
        SW  = bilread(hs3s_tile(type,root, Northing-5, Easting,CONT  ), Slat, Wlon);
        SE  = bilread(hs3s_tile(type,root, Northing-5, Easting+5,CONT), Slat, Elon);
        if NW==0 ; return; end;
        if NE==0 ; return; end;
        if SW==0 ; return; end;
        if SE==0 ; return; end;
        dem = [NW NE; SW SE];
    else % case 8; overlap near S border
        Xlon = [wwest, weast];
        N  = bilread(hs3s_tile(type,root, Northing  , Easting,CONT), Nlat, Xlon);
        S  = bilread(hs3s_tile(type,root, Northing-5, Easting,CONT), Slat, Xlon);
        if N==0; return; end;
        if S==0; return; end;
        dem = [N; S];
    end
else % cases 4,5,6
    Xlat = [wsouth, wnorth];
    if wwest < mwest % case 4; overlap near W border
        Wlon = [wwest, mwest];
        Elon = [mwest, weast];
        W  = bilread(hs3s_tile(type,root, Northing, Easting-5,CONT), Xlat, Wlon);
        E  = bilread(hs3s_tile(type,root, Northing, Easting,CONT  ), Xlat, Elon);
        if W==0; return; end;
        if E==0; return; end;
        dem = [W E];
    elseif weast > meast % case 6; overlap near E border
        Wlon = [wwest, meast];
        Elon = [meast, weast];
        W  = bilread(hs3s_tile(type,root, Northing, Easting,CONT  ), Xlat, Wlon);
        E  = bilread(hs3s_tile(type,root, Northing, Easting+5,CONT), Xlat, Elon);
        if W==0; return; end;
        if E==0; return; end;
        dem = [W E];
    else % case 5; no overlap with any border
        Xlon = [wwest, weast];
        dem  = bilread(hs3s_tile(type,root, Northing, Easting,CONT), Xlat, Xlon);
    end
end

end


% % function bilname = hs3s_tile(type,path, Northing, Easting)
% % % Create Hydrosheds 3sec file name from Northing/Easting
% % 
% % if Northing<0, NSpart = sprintf('s%02d', -Northing);
% % else           NSpart = sprintf('n%02d',  Northing);
% % end
% % 
% % if Easting<0,  WEpart = sprintf('w%03d', -Easting);
% % else           WEpart = sprintf('e%03d',  Easting);
% % end;
% % 
% % bilname = sprintf('%s\\data\\NAM\\bil\\3s\\bilcon\\%s%s_%s_bil\\%s%s_%s.bil',path, NSpart, WEpart, type,NSpart, WEpart, type);
% % end

function name = hs3s_tile(type,path, N, E,CONT)

if N<0,
    NS = 's'; N = -N;
else
    NS = 'n';
end
if E<0,
    WE = 'w'; E = -E;
else
    WE = 'e';
end;

if E<10
    if N<10
        name = sprintf('%s\\data\\%s\\bil\\3s\\bil%s\\%s0%d%s00%d_%s_bil\\%s0%d%s00%d_%s.bil',path,CONT,type,NS,N,WE,E,type,NS,N,WE,E,type);
    else
        name = sprintf('%s\\data\\%s\\bil\\3s\\bil%s\\%s%d%s00%d_%s_bil\\%s%d%s00%d_%s.bil',path,CONT,type, NS,N,WE,E,type, NS,N,WE,E,type);
    end
elseif N<10
    if E<100
        name = sprintf('%s\\data\\%s\\bil\\3s\\bil%s\\%s0%d%s0%d_%s_bil\\%s0%d%s0%d_%s.bil',path,CONT,type, NS,N,WE,E,type, NS,N,WE,E,type);
    else
        name = sprintf('%s\\data\\%s\\bil\\3s\\bil%s\\%s0%d%s%d_%s_bil\\%s0%d%s%d_%s.bil',path,CONT,type, NS,N,WE,E,type, NS,N,WE,E,type);
    end
elseif E<100;
    name = sprintf('%s\\data\\%s\\bil\\3s\\bil%s\\%s%d%s0%d_%s_bil\\%s%d%s0%d_%s.bil',path,CONT,type, NS,N,WE,E,type, NS,N,WE,E,type);
else name = sprintf('%s\\data\\%s\\bil\\3s\\bil%s\\%s%d%s%d_%s_bil\\%s%d%s%d_%s.bil',path,CONT,type, NS,N,WE,E,type, NS,N,WE,E,type);
end;

end