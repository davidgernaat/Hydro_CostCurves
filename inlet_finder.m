function inlets = inlet_finder(Z,Q,acc,fdir,WDPA_PL10,WDPA_PL20,protarea_constr,nr,nc,Zout,ZPin,deselect_mainstream,index_inlet_win,adir,flowdist,sradius)

% Z=Z_inlet_win{k};
% Q=Q_inlet_win{k};
% acc=acc_inlet_win{k};
% fdir=fdir_inlet_win{k};
% WDPA_PL10=WDPA_PL10_inlet_win{k};
% WDPA_PL20=WDPA_PL20_inlet_win{k};
% protarea_constr=1;
% nr=nr_ZiW;
% nc=nc_ZiW;
% Zout;
% ZPin=Zin(l);

defdirs

%upstream = fastfindupstream(acc,fdir,drow,dcol,Zout);
upstream = fastfindupstream_Dis(acc,fdir,flowdist,drow,dcol,Zout,sradius); %sradius search radius 50 cells = 25km
    

inlets = logical(zeros(nr,nc, 'single'));

% Selection of inlets based on settings
if protarea_constr==1
    locs = find(upstream & Z>=ZPin & Q>0.2 & WDPA_PL10==0 & WDPA_PL20==0)';
else
    locs = find(upstream & Z>=ZPin & Q>0.2)';
end

if deselect_mainstream==1
    mainstream = find_mainstream(Q,index_inlet_win,adir);
    locs(ismember(locs,mainstream))= NaN;
    locs = locs(find(~isnan(locs)));
end

for i = locs
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    % find downstream nbr
    rj = r+drow(fdir(i));
    cj = c+dcol(fdir(i));
    inlets(r,c) = Z(rj,cj) < ZPin;
end

end