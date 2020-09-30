function [r_damst,c_damst,outIdx,latOut,lonOut, grand_id] = findGrandDams(root,Rw,Q)

%% Read in GrandDB for the globe
fname = sprintf('%s\\data\\data_prep\\Grand\\grand_dams_latlonDH.csv', root);
fileID = fopen(fname);
C = textscan(fileID,'%s %s %s %s %s %s','Delimiter',';','HeaderLines',1);
% Headers in the csv: GRAND_ID	DAM_HGT_M	DAM_LEN_M	DIS_AVG_LS	LONG_DD	LAT_DD
fclose(fileID);

for j=1:numel(C)
    for i=1:numel(C{3})
        dams(i,j) = str2num(C{j}{i})';
    end
end

%% Add important dams that are missing in the Grand database

% Inga
dams(end+1,6)=   -5.531517; %lat
dams(end,5)=     13.619522; %lon
dams(end,2)=     96; %Dam height

% Bonneville
dams(end+1,6)=  45.6466; %lat
dams(end,5)=    -121.9271;%lon
dams(end,2)=    60; %dam height


%% Convert locations to r,c and select only hydroelectric dams within the basin latlon limits
for i=1:numel(dams(:,6))
    if dams(i,5) > Rw.Lonlim(1) & dams(i,5) < Rw.Lonlim(2) & dams(i,6) > Rw.Latlim(1) & dams(i,6) < Rw.Latlim(2)
        
        [r_dams(i), c_dams(i)] = setpostn(Q, Rw, dams(i,6), dams(i,5)); %Convert latitude-longitude to data grid rows and columns
        grand_id(i)=dams(i,1);
    else
        r_dams(i)=NaN; c_dams(i)=NaN; grand_id(i)=NaN;
    end
end

%% remove NaNs first time
r_dams=r_dams((~isnan(r_dams)));
c_dams=c_dams((~isnan(c_dams)));
grand_id=grand_id((~isnan(grand_id)));

%% Remove locations outside Basin
for i=1:numel(r_dams)
    if r_dams(i)==0; continue; end;
    if Q(r_dams(i),c_dams(i))==0; r_dams(i)=NaN; c_dams(i)=NaN; grand_id(i)=NaN; end
end

%% remove NaNs 2nd time
r_dams=r_dams((~isnan(r_dams)));
c_dams=c_dams((~isnan(c_dams)));
grand_id=grand_id((~isnan(grand_id)));

%%
if isempty(r_dams)==1;
    r_damst=-99;c_damst=-99;outIdx=-99;latOut=-99;lonOut=-99;
    disp('findGrandDams: no dams in current basin')
    return
end

%%
%     figure(1);clf;imagesc(log(Q)); axis image; colormap(flipud(gray));
%     hold on;
%     pDams=plot(c_dams,r_dams,'b.','markersize',15); %Existing dams
%     hold off;

%% Correct GrandDB locations to cells with high Q
%
disp('findGrandDams: Correcting locations')
for i=1:numel(r_dams)
    Qdams_preCorr(i) = Q(r_dams(i),c_dams(i));
end

%%
counter=0;
for i=1:numel(r_dams)
    win=1;
    firstrow = r_dams(i)-win;
    lastrow  = r_dams(i)+win;
    firstcol = c_dams(i)-win;
    lastcol  = c_dams(i)+win;
    Qwin = Q(firstrow:lastrow,firstcol:lastcol); %Select small Q window
    
    Qfocus = Q(r_dams(i),c_dams(i));
    Qmax = max(Qwin(:));
    Qratio = Qmax/Qfocus;
    
    if Qratio>2;
        counter=counter+1;
        %fprintf('Correction #%d\n',counter);
        [r_max,c_max]=find(Qwin==max(Qwin(:))); %Find max
        
        %Transform r,c to highest point
        r_trans = r_max(1)-win-1; %In case multiple max, first option
        c_trans= c_max(1)-win-1;
        
        %Transform to large grid
        r_damst(i) = r_dams(i)+r_trans;
        c_damst(i) = c_dams(i)+c_trans;
        
        %check
        %figure(1);clf;imagesc(Qwin);axis image;colormap(flipud(gray)); hold on
        %plot(win+1,win+1,'r.','markersize',20);
        %plot(win+1+c_trans,win+1+r_trans,'b.','markersize',20);
        %hold off
    else
        r_damst(i) = r_dams(i);
        c_damst(i) = c_dams(i);
    end
end

%%
for k=1:numel(r_damst)
    
    if isnan(r_damst(k))==1; continue; end;
    
    outIdx(k) = sub2ind(size(Q),r_damst(k),c_damst(k));
    
end

%% New corrected coordinates
for k=1:numel(r_dams)
    [latOut(k), lonOut(k)] = setltln(Q, Rw, r_damst(k), c_damst(k)); %Coordinates Grand dams
end

end