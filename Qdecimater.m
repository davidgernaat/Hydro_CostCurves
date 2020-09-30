function QQ = Qdecimater(rinMin,cinMin,inlet_win,Q,fdir,acc,inletOpt)

% rinMin=84;
% cinMin=81;
% inlet_win=75;
% Q = Q_inlet_win{k};
% fdir=fdir_inlet_win{k};
% acc=acc_inlet_win{k};
% inletOpt=aPinlet_windowMin(k);

stream = logical(zeros(size(Q)));
QQ = Q;

r = rinMin;
c = cinMin;

if isempty(rinMin)==1; return; end;
if isempty(cinMin)==1; return; end;

defdirs

while 1
    [r c];
    stream(r,c) = true;
    d = fdir(r,c);
    r = r+drow(d);
    c = c+dcol(d);
    if r==inlet_win+1 && c==inlet_win+1, break; end; %Stop at outlet
end;

dQ = QQ(rinMin,cinMin);
QQ(stream) = max(0,QQ(stream) - dQ);

% Also remove Q upstream of inlet
upstream = fastfindupstream_lim(acc,fdir,drow,dcol,inletOpt);
QQ(upstream) = 0;

end