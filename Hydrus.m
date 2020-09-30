tic
%% Backbone Hydropower model - Basic settings
%  Convert to numbers
if isdeployed
    nbasin_start = str2num(nbasin_start);
    nbasin_end = str2num(nbasin_end);
    protarea_constr = str2num(protarea_constr);
    navi_constr = str2num(navi_constr);
    ExistingDams_constraint = str2num(ExistingDams_constraint);
    BeforeFirstDam_constraint = str2num(BeforeFirstDam_constraint);
    slackflow_constraint = str2num(slackflow_constraint);
    do = str2num(do);
    betaL = str2num(betaL);
    waterconsumption = str2num(waterconsumption);
end
if ~isdeployed
    clc
    clear all
    root = pwd;%'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\Hydropower\Model\Sanita_model_package';
    root_bil = pwd;% 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\Hydropower\Model\Sanita_model_package'; %For grid calculations
    continent_in = 'ASIA';
    nbasin_start = 5;
    nbasin_end = 5;
    version='Remain';
    
    %Settings
    protarea_constr=1;           %No dams in protected areas
    navi_constr=0;               %No dams in navigable rivers (see depth_cutoff-setting)
    ExistingDams_constraint=1;   %No dams on locations in existing hydro lakes
    BeforeFirstDam_constraint=1; %No dams on mainstream locations before the first GrandDam
    slackflow_constraint=0;      %Slackflow of x% to ensure natural flow of the river (see slackflow -setting)
    do=5000;                       %Distance between the outlets (spacing), in terms of number of cells. 1 cell is ~500m. Standard:50 (25km); for development:2000(1000km)
    betaL=0;                     %Deselection based on Lake size. Cost per power is always on. This setting should be on only in eco scenario
    waterconsumption=0;          %Discharge with (1) or without water consumption (0)
    nbasin = nbasin_start;
end

%% Log all outputs
runname=sprintf('do%d_woWC_TEST00_Original',do);
matpath = fullfile(root, sprintf('output\\%s\\%s\\%s', version, continent_in,runname));
if ~isfolder(matpath); mkdir(matpath); end

logf = fullfile(matpath, sprintf('Hydrus_%s_%s.log',string(datetime('now','Format','dd-MMM-yyyy')),runname));
diary(logf)

%enable profiler to track runtime and dependencies
saveProfile=1;
if saveProfile; profile on; end

%% Start run
mystartstatement='\n HYDRUS00 START TEST RUN for do=%d:  %s\n';
fprintf(mystartstatement,do,datetime('now'));

for nbasin = nbasin_start:nbasin_end
    clearvars -except root root_bil continent_in nbasin nbasin_start nbasin_end version protarea_constr navi_constr ExistingDams_constraint ...
        BeforeFirstDam_constraint slackflow_constraint do betaL waterconsumption matpath saveProfile
    
    disp('Start')
    
    %     try %to prevent loop from crashing despite error
    
    fprintf('\nSETTINGS\n\n');
    fprintf('Continent: %s\n', continent_in);
    fprintf('Basin: %d of %d\n\n', nbasin, nbasin_end);
    fprintf('Scenario: %s\n', version);
    fprintf('Protected areas constraint: %d\n', protarea_constr);
    fprintf('Navigational constraint: %d\n', navi_constr);
    fprintf('Existing dams constraint: %d\n', ExistingDams_constraint);
    fprintf('No dams before first dam: %d\n', BeforeFirstDam_constraint);
    fprintf('Slack flow constraint: %d\n', slackflow_constraint);
    fprintf('Spacing distance: %d\n', do);
    fprintf('Deselection based on lake size: %d\n', betaL);
    fprintf('Discharge corrected for water consumption: %d\n\n', waterconsumption);
    
    % Reading ancillary scripts and data
    disp('Reading ancillary scripts and data');
    % add functions
    if ~isdeployed
        addpathname = sprintf('%s\\functions', root);
        addpath(addpathname)
    end
    
    % add routing script
    defdirs
    
     % First data file for current basins
    fname = sprintf('%s\\data\\%s\\Basin\\Basin_%d.mat', root, continent_in, nbasin);
    load(fname);
    [nrw,ncw] = size(acc);
    
    % Second all basins in continent
    fname = sprintf('%s\\data\\%s\\basins.mat', root_bil, continent_in);
    load(fname);
    [nr,nc] = size(basin);
    
    % Third reading GDPpc ppp
    fname = sprintf('%s\\data\\data_prep\\GDPpc\\GDPpc2010IsoCode.csv', root_bil);
    fileID = fopen(fname);
    C = textscan(fileID,'%s %s %s %s %s %s %s','Delimiter',';','HeaderLines',1);
    fclose(fileID);
    for i=1:numel(C{3})
        ISOGDP(i,1) = str2num(C{3}{i})'; %ISO number
        ISOGDP(i,2) = str2num(C{4}{i})'; %GDP pc
        ISOGDP(i,3) = str2num(C{5}{i})'; %GDP MER
    end
    
    %% Create Indus basin window as done in basin_selector.m for Georef
    Z=single(Z);
    desbasin = basin==nbasin; %Desired Basin, 4 is Columbia
    
    % get columns of basin for lon lat coordinates
    cols = find(sum(desbasin)>0); % columns
    rows = find(sum(desbasin')>0); % rows
    
    % Window around basin, comes from basin_selector.m
    % Important for geolocation
    georefwin=minwin;
    
    % Large window size to circumvent problems of winsize down the model
    % Careful! Needs to be consistent with winsize basine_selection.m
    firstcol = max(1,cols(1) - georefwin);
    lastcol = min(nc,cols(end) + georefwin);
    firstrow = max(1,rows(1) - georefwin);
    lastrow = min(nr,rows(end) + georefwin);
    
    % Additional windows
    [~,orderacc] = sort(acc(:),'descend');  %For fastfindupstream
    inlets = zeros(nrw,ncw, 'int8');        %For find_inlets
    
    % Set Geo reference coordinates and settings
    georef % Creates the Rw georef structure used throughout to go from ltln to rc and viceversa
    
    returnID=0;
    %% Set output file locations
    %matpath = fullfile(root, sprintf('output\\%s\\%s\\do%d_woWC_TEST00', version, continent_in, do));
    matfile = fullfile(matpath, sprintf('Basin%d_output_do%d.mat', nbasin,do));
    matfileCOEPOT = fullfile(matpath, sprintf('COEPOT_b%d_do%d.mat', nbasin,do));
    if ~isfolder(matpath); mkdir(matpath); end
    
    %% Get GrandDatabase location
    [r_damst,c_damst,GrandIdx,Grandlat,Grandlon] = findGrandDams(root_bil,Rw,Q);
    
    %% Run Settings
    
    % RD = River power systems (River dams)
    % DP = Reservoir systems (Damp-Pipe systems)
    % P  = Diversional canal systems (Pipe systems)
    
    % k = loop for outlet locs
    % l = loop for inlet elevation levels
    % m = loop for inlet per elevation level
    % ndl = loop for Q decreaser
    
    nd=4;                        %Number of Q-decreaser loops
    dowin=50;                    %Winsize for DP and P calculations.
    sradius=40;                  %Search radius for DP and P systems
    ni=6;                        %Number of inlet elevation levels per outlet loc
    cost_lim=0.5;                %Cost limit $/kWh
    outletdeselect=2;            %Number of deselected sites upstream of outlet
    rivermouth_inland=400;       %Number of cells inland of basin outlet deselected (200km) (rivermouth_constraint)
    depth_cutoff=4;              %Water depth cutoff (m)
    bighydro_cutoff=50;          %Max Q based on small-scale Hydropower Veileder 10MW document p20
    MiniDamMax=30;               %Max Dam Height based on Q cutoff. This is consistent with big Hydro Veileder p62
    slackflow=0.3;               %Percentage of water diverted around the dam to ensure natural ecological flow of the river
    seismicpikcost=1;            %Including seismic cost yes(1) or no(0)
    
    % Constraints. 0=constraint off / 1= constraint on
    runDP=0;                     %Run dam-pipe-resevoir systems 0=no/1=yes
    runP=1;                      %Run pipe systems 0=no/1=yes
    runRD=1;                     %Run river dam system 0=no/1=yes
    runMini=0;                   %Run mini dam systems 0=no/1=yes
    
    cost_constr=1;               %Enabling cost constraint (see cost_lim-setting)
    rivermouth_constr=1;         %No dams in rivermouth (see outletdeselect-setting)
    mainstream_constr=0;         %Deselect dam locations on basin mainstream
    deselect_mainstream=0;       %Deselect DP and P locations that are on the mainstream of RD power station.
    MiniHydro_deselect=0;        %Deselect dam locations on small streams (MiniHydro locations) (see bighydro_cutoff-setting)
    MiniHydro_select=0;          %Select dam locations on small streams (MiniHydro locations) (see bighydro_cutoff-setting)
    MiniHydro_special=0;         %Special costmodel activation for MiniHydro locations (see bighydro_cutoff-setting and MiniDamMax-setting)
    mangrove_constr=0;           %No dams in mangroves
    seismicnoaa_constr=0;        %No dams in earth quake areas NOAA
    discost_set=1;               %Taking into account distance to load cost. 0=no/1=yes
    landval_set=1;               %Taking into account Land value cost. 0=no/1=yes
    ASYMRD=0;                    %Asymmetric dams RD-systems 0=no/1=yes
    ASYMDP=0;                    %Asymmetric dams DP-systems 0=no/1=yes
    betaC=0;                     %Deselection cost weight
    betaP=0;                     %Deselection power weight
    betaCP=1;                    %Deselection cost per power weight
    betaA=0;                     %Deselection flow acc (as indication of river lenght)
    
    %% Find outlets and create windows for Qdecreaser
    disp('Assigning powerstation locations');
    
    %%%%find_outlets_wins.m START
    %% Skip if max flowdist is < do (distance between outlets)
    
    if max(flowdist(:)) < do;
        disp('Flowdist is less than distance between outlets');
        save(matfile,'-v7.3','do');
        returnID=1;
        return
    end
    
    %% Find outlets and create windows
    
    [ro,co] = find_outlets_dis(flowdist,acc,do);
    outIdx = sub2ind(size(acc),ro,co);
    
    %Make variables
    for k=1:numel(ro)
        Accoutlets(k)=0;
        PopDisplacedOpt(k)=0;
        Quakerate(k)=0;
        COETotRD(k)=0;
        RDPnet(k)=0;
        Qoutlets_design(k)=0;
        Qoutlets_design_LF(k)=0;
        RDRegion_id(k)=0;
        RDCountry_id(k)=0;
        latOut(k)=0;
        lonOut(k)=0;
        OptDH(k)=0;
        OptDL(k)=0;
        OptPop(k)=0;
        OptLV(k)=0;
        RDDepth(k)=0;
        OptInv(k)=0;
        RDP(k)=0;
        RDlakeSurface{k}=[];
        RDVolumeLake{k}=[];
        RDVolumeLake15s(k)=0;
        RDSurfaceLake15s(k)=0;
        Zoutlets(k)=0;
        DisOutlet(k)=0;
        OptSpecCap(k)=0;
        PopDisplacedOpt(k)=0;
        Zoutlets(k)=0;
        Qoutlets(k)=0;
        Qoutlets_design(k)=0;
        Qoutlets_design_LF(k)=0;
        RDRegion_id(k)=0;
        RDCountry_id(k)=0;
        RDDepth(k)=0;
        Accoutlets(k)=0;
        
    end
    
    % figure(1);clf;imagesc(log(Q));colormap(flipud(gray));axis image; hold on
    % plot(co,ro,'r.','markersize',15); hold off
    % sum(isnan(ro))
    % figure(2);ax2=subplot(1,3,2);imagesc(log(Q));colormap(flipud(gray));axis image; hold on
    % plot(co,ro,'r.','markersize',15); hold off
    %% Q-maps with ot without water consumption
    
    if waterconsumption==1
        Q=Qwc;
        Qdesign=Qdesign_wc;
        Qdesign_mean=Qdesign_mean_wc;
        Qdesign_LF=Qdesign_LF_wc;
    end
    
    %% Lat lon coordinates outlets
    
    for  k = 1:numel(ro)
        [latOut(k), lonOut(k)] = setltln(acc, Rw, ro(k), co(k)); %Coordinates outlets
        %     latOut(k) = Rlatw(ro(k))-((Rlatw(2)-Rlatw(1))/2);
        %     lonOut(k) = Rlonw(ro(k))-((Rlonw(2)-Rlonw(1))/2);
        Depth_Rivermouth(k) = D(ro(k),co(k)); % Depth of locations from outlet
        flowdist_rivermouth(k) = flowdist(ro(k),co(k)); % Distance from outlet
        Qoutlets(k) = Q(ro(k),co(k));
        outlets(k) = sub2ind(size(Z),ro(k),co(k));
        Accoutlets(k) = acc(ro(k),co(k));
    end
    
    %% Deselect outlet and - No dam in river mouth
    
    if rivermouth_constr==1
        if numel(ro)<outletdeselect;
            disp('To few ro and co');
            save(matfile,'-v7.3','do');
            returnID=1;
            return
        end
    end
    
    %% Rivermouth constraint
    if rivermouth_constr==1
        
        [~,sortQidx]=sort(Accoutlets,'descend');
        
        %First two for sure
        ro(sortQidx(1:outletdeselect))=NaN;
        co(sortQidx(1:outletdeselect))=NaN;
        %     ro(sortQidx(50:(end)))=NaN;
        %     co(sortQidx(50:(end)))=NaN;
        
        %This method causes additional deselection after lakes in the river
        %system, especially in Congo and Amazon
        %Second deselect dams on mainstream xkm inland and on rivers with a certain depth
        for k=1:numel(ro)
            if flowdist_rivermouth(k) < rivermouth_inland & Depth_Rivermouth(k) > depth_cutoff
                ro(k)=NaN;
                co(k)=NaN;
            end
        end
        
        %Other 2 maybe depending on river depth
        %     if strcmp('SAM',continent_in)
        %         if nbasin==1
        %             if numel(sortQidx) < 8; RDi = numel(sortQidx); else RDi=33; end
        %             for i=3:RDi
        %                 if Depth_Rivermouth(sortQidx(i)) > depth_cutoff
        %                     ro(sortQidx(i))=NaN;
        %                     co(sortQidx(i))=NaN;
        %                 end
        %             end
        %         end
        %     end
    end
    % sum(isnan(ro))
    
    % figure(1); hold on
    % plot(co,ro,'b.','markersize',15); hold off
    % figure(2); subplot(1,3,2); hold on
    % plot(co,ro,'b.','markersize',15); hold off
    % fd=flowdist;fd(Q>100)=0;cmap=jet(256);cmap(1,:)=[1 1 1];
    % figure(2);ax3=subplot(1,3,3); imagesc(fd);colormap(cmap);axis image;
    
    %% Deselect locations in WDPA regions
    
    if protarea_constr==1
        for k= 1:numel(ro)
            if isnan(ro(k))==1; continue; end;
            if WDPA_PL10(ro(k),co(k))== 10 || WDPA_PL20(ro(k),co(k))==20
                ro(k)=NaN;
                co(k)=NaN;
            else
            end
        end
    else
    end
    
    %% Navigation constraint for navigability (4m)
    
    if navi_constr==1
        for k= 1:numel(ro)
            if isnan(ro(k))==1; continue; end;
            if D(ro(k),co(k))> depth_cutoff
                ro(k)=NaN;
                co(k)=NaN;
            else
            end
        end
    else
    end
    % [rd, cd] = setpostn(Q, Rw, -3.128930, -60.027996);
    % figure(2); hold on
    % plot(cd,rd,'g.','markersize',25);
    % plot(co,ro,'b.','markersize',15); hold off
    % sum(isnan(ro))
    %% Mangrove constraint
    
    if mangrove_constr==1
        for k= 1:numel(ro)
            if isnan(ro(k))==1; continue; end;
            if MC(ro(k),co(k))==1
                ro(k)=NaN;
                co(k)=NaN;
            end
        end
    end
    
    %% Mainstream-constraint - No dams on basin mainstream
    
    if mainstream_constr==1
        [~,MainIdx] = max(Q(:));
        mainstrm = find_mainstream(Q,MainIdx,adir);
        outlets(ismember(outlets,mainstrm))= NaN;
        ro(find(isnan(outlets)))=NaN;
        co(find(isnan(outlets)))=NaN;
    end
    
    %% seismicnoaa-constraint - No dams in earthquake areas
    
    if seismicnoaa_constr==1
        for k= 1:numel(ro)
            if isnan(ro(k))==1; continue; end;
            if seismicnoaa(ro(k),co(k))==1
                ro(k)=NaN;
                co(k)=NaN;
            end
        end
    end
    
    %% MiniHydro deselection
    
    if MiniHydro_deselect==1
        for k= 1:numel(ro)
            if isnan(ro(k))==1; continue; end;
            if Q(ro(k),co(k))< bighydro_cutoff
                ro(k)=NaN;
                co(k)=NaN;
            end
        end
    end
    
    %% MiniHydro Selection
    
    if MiniHydro_select==1
        for k= 1:numel(ro)
            if isnan(ro(k))==1; continue; end;
            if Q(ro(k),co(k))> bighydro_cutoff
                ro(k)=NaN;
                co(k)=NaN;
            end
        end
    end
    
    %% Existing dams de-selection
    
    if ExistingDams_constraint==1
        for k= 1:numel(ro)
            if isnan(ro(k))==1; continue; end;
            if NoDamsLand(ro(k),co(k))== 1
                ro(k)=NaN;
                co(k)=NaN;
            end
        end
    end
    
    %% No dams before first GrandDam on mainstream
    
    if BeforeFirstDam_constraint==1
        if r_damst(1)~=99;
            NoBeforeFirstMainstrDam = FirstMainstrmDam(Q,adir,acc,r_damst,c_damst,GrandIdx);
            outlets(ismember(outlets,NoBeforeFirstMainstrDam))= NaN;
            ro(find(isnan(outlets)))=NaN;
            co(find(isnan(outlets)))=NaN;
        end
    end
    
    %% Deselect dam location with zero discharge value (can happen due to water consumption)
    
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if Q(ro(k),co(k))== 0
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
    
    %% if all NaN break loop
    
    if sum(isnan(ro))==numel(ro);
        disp('First all NaN break');
        save(matfile,'-v7.3','do');
        returnID=1;
        return
    end
    
    %% Transform
    for k=1:numel(ro)
        
        if isnan(ro(k))==1; continue; end;
        
        Zoutlets(k) = Z(ro(k),co(k));
        Qoutlets(k) = Q(ro(k),co(k));
        Qoutlets_design(k) = Qdesign(ro(k),co(k));
        Qoutlets_design_LF(k) = Qdesign_LF(ro(k),co(k));
        RDRegion_id(k) = Regions(ro(k),co(k));
        RDCountry_id(k) = Countries(ro(k),co(k));
        RDDepth(k) = D(ro(k),co(k));
        DisOutlet(k) = Dis(ro(k),co(k)); % Distance to powerline (km)
    end
    
    %% Slackflow for ecological reasons
    if slackflow_constraint==1
        Qoutlets_design = Qoutlets_design *( 1-slackflow);
    end
    
    %% Control for NaNs in Dis map (slightly inconsistent with other maps because different sea map was used)
    for k=1:numel(DisOutlet)
        if isnan(DisOutlet(k))==1; ro(k)=NaN; co(k)=NaN; end;
    end
    
    % figure(2); hold on
    % plot(co,ro,'b.','markersize',15); hold off
    
    %% Seismic PIK
    % 5 percent cost increase at seismic hazard > 4
    
    if seismicpikcost==1
        for k=1:numel(ro)
            if isnan(ro(k))==1; Quakerate(k)=NaN; continue; end;
            
            if seismicpik(ro(k),co(k)) > 4
                Quakerate(k) = 0.05;
            end
        end
    else
        for k=1:numel(ro)
            Quakerate(k) = 0;
        end
    end
    
    %% Report deselection
    no=numel(outlets);
    noNAN=sum(isnan(ro));
    fprintf('Number of outlets: %d.\n', no);
    fprintf('Number of deselected: %d.\n', noNAN);
    
    %% Create windows for Qdecreaser
    
    for k = 1:numel(outlets)
        
        if isnan(ro(k))==1; continue; end;
        
        inlet_win = dowin;
        ZiWinfr(k) = ro(k)-inlet_win;
        ZiWinlr(k) = ro(k)+inlet_win;
        ZiWinfc(k) = co(k)-inlet_win;
        ZiWinlc(k) = co(k)+inlet_win;
        
        if ZiWinfr(k)<1 || ZiWinlr(k)>nrw || ZiWinfc(k)<1 || ZiWinlc(k)>ncw
            ro(k)=NaN;
            co(k)=NaN;
        else
        end
        
        if isnan(ro(k))==1; continue; end;
        
        
        fdir_inlet_win{k}      = fdir(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        adir_inlet_win{k}      = adir(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        Q_inlet_win{k}         = Q(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        acc_inlet_win{k}       = acc(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        Z_inlet_win{k}         = Z(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        Pop_inlet_win{k}       = Pop(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        WDPA_PL10_inlet_win{k} = WDPA_PL10(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        WDPA_PL20_inlet_win{k} = WDPA_PL20(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        LandValue_inlet_win{k} = LandValue(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        Dis_inlet_win{k}       = Dis(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        flowdist_inlet_win{k}  = flowdist(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
        
        
    end
    
    WinProb=sum(isnan(ro))- noNAN;
    fprintf('Number of window problems: %d.\n', WinProb);
    
    index_inlet_win = sub2ind([(2*inlet_win+1) (2*inlet_win+1)],inlet_win+1,inlet_win+1);
    
    %% if all NaN break loop
    
    if sum(isnan(ro))==numel(ro);
        disp('all NaN brak loop');
        save(matfile,'-v7.3','do');
        save(matfileCOEPOT,'-v7.3','do');
        returnID=1;
        return
    end
    
    % figure(2); hold on
    % plot(co,ro,'b.','markersize',15); hold off
    
    %%
    % disp('Saving Depth')
    % matpath = fullfile(root, sprintf('output\\%s\\%d', continent_in, nbasin));
    % matfile = fullfile(root, sprintf('output\\%s\\%d\\Basin%d_Depth.mat', continent_in, nbasin, nbasin));
    % if ~isdir(matpath)
    %     mkdir(matpath);
    % end
    %
    % save(matfile,'-v7.3','RDDepth');
    %
    % continue
    
    % disp('Test');
    % save(matfile,'-v7.3','do','ro','co','ro_org','co_org');
    % returnID=1;
    % return
    %%%%find_outlets_wins.m END
    
    %
    %figure(1);clf;imagesc(log(Q));axis image;colormap(flipud(gray));hold on
    %figure(1);clf;imagesc(seismicpik);axis image;colormap(parula);hold on
    %plot(c_damst,r_damst,'.r','markersize',20)
    %plot(co,ro,'.b','markersize',20)
    %hold off
    
    % catch error
    %    if returnID; continue; end
    
    %% RUN-OF-RIVER PIPE-SYSTEMS, NO DAM, LOADFACTOR <100%
    if runP==1
        for ndl = 1:nd
            
            fprintf('Pipe Systems Q decreaser: %d of %d.\n', ndl, nd);
            
            %%%% PipeSystems.m START
            %% Pipe-systems procedure
            
            %% Find pipe inlets
            disp('Find Pinlets');
            
            %%%%find_Pinlets.m START
            %% Find Pinlet locations
            
            for k = 1:numel(outlets);
                
                if isnan(ro(k))==1; Pinlet_win{k}=[]; continue; end;
                
                %Data windows are created in find_outlet_wins
                
                [nr_ZiW,nc_ZiW] = size(Z_inlet_win{k});
                
                Z_inlet_win_temp=0;
                rr1=inlet_win+1-sradius;
                rr2=inlet_win+1+sradius;
                cc1=inlet_win+1-sradius;
                cc2=inlet_win+1+sradius;
                Z_inlet_win_temp = Z_inlet_win{k}(rr1:rr2,cc1:cc2);
                
                Zmin = Z_inlet_win{k}(inlet_win+1,inlet_win+1);
                Zmax = max(Z_inlet_win_temp(:));
                Zupdown = Zmax-Zmin;
                ZPin = linspace(Zmin+(0.1*Zupdown),Zmax-(0.1*Zupdown),ni);
                
                Zout=  sub2ind(size(Z_inlet_win{k}),inlet_win+1,inlet_win+1); %index of center dam outlet location
                
                %figure(2);clf;imagesc(Z_inlet_win{k});axis image; colormap(flipud(gray));hold on
                %plot(inlet_win+1,inlet_win+1,'.r','markersize',20); hold off
                
                for l = 1:numel(ZPin)
                    
                    %Pinletsw = inlet_finder(Z_inlet_win{k},Q_inlet_win{k},acc_inlet_win{k},fdir_inlet_win{k},WDPA_PL10_inlet_win{k},WDPA_PL20_inlet_win{k},protarea_constr,nr_ZiW,nc_ZiW,Zout,ZPin(l),deselect_mainstream,index_inlet_win,adir_inlet_win{k});
                    Pinletsw = inlet_finder(Z_inlet_win{k},Q_inlet_win{k},acc_inlet_win{k},fdir_inlet_win{k},WDPA_PL10_inlet_win{k},WDPA_PL20_inlet_win{k},protarea_constr,nr_ZiW,nc_ZiW,Zout,ZPin(l),deselect_mainstream,index_inlet_win,adir_inlet_win{k},flowdist_inlet_win{k},sradius);
                    
                    Pinlet_win{k}{l} = find(Pinletsw);
                    %convert inlet maps to arrays of rows/columns
                    rPin_w{k}{l} = rem(Pinlet_win{k}{l}-1, nr_ZiW)+1;
                    cPin_w{k}{l} = fix((Pinlet_win{k}{l}-1)/nr_ZiW)+1;
                    % Converting r,c to Columbia catchment r,c
                    rPin{k}{l} = ro(k)+(rPin_w{k}{l}-inlet_win-1);
                    cPin{k}{l} = co(k)+(cPin_w{k}{l}-inlet_win-1);
                    Pinlet{k}{l} = sub2ind(size(Z),rPin{k}{l},cPin{k}{l});
                    ZPinlet{k}{l} = Z(Pinlet{k}{l});
                    QPinlet{k}{l} = Q(Pinlet{k}{l});
                    accPinlet{k}{l} = acc(Pinlet{k}{l});
                    QDesignPinlet{k}{l} = Qdesign(Pinlet{k}{l});
                    QDesignLFPinlet{k}{l} = Qdesign_LF(Pinlet{k}{l});
                    QDesignMeanPinlet{k}{l} = Qdesign_mean(Pinlet{k}{l});
                    PRegion_id{k}{l} = Regions(rPin{k}{l},cPin{k}{l});
                    PCountry_id{k}{l} = Countries(rPin{k}{l},cPin{k}{l});
                    
                    
                    if slackflow_constraint==1
                        QDesignPinlet{k}{l} = QDesignPinlet{k}{l} *( 1-slackflow);
                    end
                    
                    %Check
                    %figure(3);clf;imagesc(log(Q_inlet_win{3}));axis image; colormap(flipud(gray));hold on
                    %plot(inlet_win+1, inlet_win+1,'b.','markersize',24);
                    %plot(cPin_w{3}{1}, rPin_w{3}{1},'r.','markersize',24);
                    %hold off
                    
                    %figure(3);clf;imagesc(log(Q));axis image; colormap(flipud(gray));hold on
                    %plot(cPin{3}{1}, rPin{3}{1},'r.','markersize',24);
                    %hold off
                    
                end
            end
            
            %% Variable creator
            
            
            for k = 1:numel(outlets)
                
                COEP{k}=[];
                PL{k}=[];
                latP{k}=[];
                lonP{k}=[];
                PPnet{k}=[];
                PPipeDia{k}=[];
                OptInvP{k}=[];
                PP{k}=[];
                bMinElevP(k)=0;
                aCOEPmin(k)=0;
                aPPnetmin(k)=0;
                aCOEPmin(k)=NaN;
                aPPnetmin(k)=NaN;
                CostElementsPMin{k}=[];
                dfQPmin(k)=NaN;
                accPmin(k)=NaN;
                dfPLmin(k)=NaN;
                HeadPmin(k)=NaN;
                HeadraceLPmin(k)=NaN;
                aPInletMin(k)=NaN;
                aPinlet_windowMin(k)=NaN;
                lat_Pin_min(k)=NaN;
                lon_Pin_min(k)=NaN;
                rPinMin(k)=NaN;
                cPinMin(k)=NaN;
                QDesignPinletMin(k) =NaN;
                QDesignLFPinletMin(k) =NaN;
                QDesignMeanPinletMin(k) = NaN;
                ZPinletMin(k) =NaN;
                nPipePMin(k) =NaN;
                OptInvPMin(k) =NaN;
                PPMin(k) =NaN;
                DeselectedSites_unSortRD{k}=[];
                DeselectedGrandSites_unSortRD{k}=[];
                MiniFlag(k) = 0;
                OptInv(k) =NaN;
                OptSpecCap(k) =NaN;
                RDP(k) = NaN;
                RDVolumeLake15s(k) = NaN;
                RDSurfaceLake15s(k) = NaN;
                OptPop(k) = NaN;
                OptLV(k) = NaN;
                OptSpecCapP{k} = [];
                OpthfP{k} = [];
                OptDP{k} = [];
                OptSpecCapPMin(k) = NaN;
                OpthfPMin(k) = NaN;
                PSurfaceLake15s(k) = 1e-3;
                OptDPMin(k) = NaN;
                
                if numel(Pinlet_win{k})==0; continue; end;
                
                for l = 1:ni
                    COEP{k}{l}=[];
                    PL{k}{l}=[];
                    latP{k}{l}=[];
                    lonP{k}{l}=[];
                    PPnet{k}{l}=[];
                    PPipeDia{k}{l}=[];
                    OptInvP{k}{l}=[];
                    PP{k}{l}=[];
                    COEPminElev{k}(l)=0;
                    bMinP{k}(l)=0;
                    OptSpecCapP{k}{l}=[];
                    OpthfP{k}{l}=[];
                    OptDP{k}{l}=[];
                    
                    for m = 1:numel(Pinlet_win{k}{l})
                        
                        COEP{k}{l}(m)=NaN;
                        PL{k}{l}(m)=NaN;
                        latP{k}{l}(m)=NaN;
                        lonP{k}{l}(m)=NaN;
                        PPnet{k}{l}(m)=NaN;
                        PPipeDia{k}{l}(m)=NaN;
                        OptInvP{k}{l}(m)=NaN;
                        PP{k}{l}(m)=NaN;
                        OptSpecCapP{k}{l}(m)=NaN;
                        OpthfP{k}{l}(m)=NaN;
                        OptDP{k}{l}(m)=NaN;
                    end
                end
            end
            
            for k=1:nd
                for l = 1:numel(outlets)
                    for m = 1:ni
                        DeselectedSites_unSortDP{k}{l}{m}=[];
                    end
                end
            end
            
            DeselectedSites = 0;
            DeselectedSites_final = {};
            
            %%%%find_Pinlets.m END
            
            %% Determine lat lon coordinates of outlets and inlets
            disp('Calculate pipelenght and lat lon');
            
            %%%%Ppipelength.m START
            %% Calculate pipeplenght
            % Based on columns rows and elevation diff
            
            for  k = 1:numel(outlets)
                
                if isnan(ro(k))==1; continue; end;
                
                for l = 1:ni
                    for m = 1:numel(Pinlet_win{k}{l})
                        
                        PL{k}{l}(m) = Lenght_finder(acc_inlet_win{k},Zout,Pinlet_win{k}{l}(m),latOut(k),Zoutlets(k),ZPinlet{k}{l}(m));
                        
                        [latP{k}{l}(m), lonP{k}{l}(m)] = setltln(acc, Rw, rPin{k}{l}(m), cPin{k}{l}(m)); %Coordinates inlets
                        
                    end
                end
            end
            %%%%Ppipelength.m END
            
            %% Running cost model
            disp('Running pipe cost model');
            
            %%%%Pcostmodel.m START
            %% Pipe cost model
            
            for  k = 1:no %outlets
                if isnan(ro(k))==1; continue; end;
                for l = 1:ni %elevations
                    for m = 1:numel(Pinlet_win{k}{l}) % inlets
                        
                        if isnan(lonP{k}{l}(m))==1;
                            COEP{k}{l}(m)= NaN;
                            PPnet{k}{l}(m) = NaN;
                            CostElementsP{k}{l}{m} = NaN;
                            nPipeP{k}{l}(m)=NaN;
                            continue; end;
                        
                        [COEP{k}{l}(m) PPnet{k}{l}(m) CostElementsP{k}{l}{m} nPipeP{k}{l}(m) OptInvP{k}{l}(m) PP{k}{l}(m) OptSpecCapP{k}{l}(m) OpthfP{k}{l}(m) OptDP{k}{l}(m)] = costmodel_pipesys(Zoutlets(k),ZPinlet{k}{l}(m),PL{k}{l}(m),QDesignPinlet{k}{l}(m),QDesignMeanPinlet{k}{l}(m),QDesignLFPinlet{k}{l}(m),DisOutlet(k),cost_constr,cost_lim,Quakerate(k));
                        
                        if COEP{k}{l}(m)==0; COEP{k}{l}(m)=NaN;end;
                        if PPnet{k}{l}(m)==0; PPnet{k}{l}(m)=NaN;end;
                        PPnet{k}{l}(m) = PPnet{k}{l}(m) *1e-6; % kWh to GWh
                        
                    end
                end
            end
            
            %% Select Cheapest
            % First find lowest COE per elevation level
            for k=1:no
                
                if isnan(ro(k))==1; continue; end;
                
                for l=1:ni
                    if isempty(COEP{k}{l})==1; continue; end;
                    [COEPminElev{k}(l),bMinP{k}(l)] = min(COEP{k}{l}); %Lowest index per elevation level
                end
            end
            
            for k=1:no
                
                if isnan(ro(k))==1; continue; end;
                
                bMinP{k}(find(bMinP{k}==0))=NaN;
                COEPminElev{k}(find(COEPminElev{k}==0))=NaN;
            end
            
            % Secondly find lowest COE among elevation levels
            for k=1:no
                
                if isnan(ro(k))==1; continue; end;
                if isempty(COEPminElev{k})==1; continue; end;
                
                
                [~,bMinElevP(k)] = min(COEPminElev{k}); %Lowest index of all elevation levels
            end
            
            %% Transfer min info into actual data
            clear a1 a2
            for k= 1:no
                
                if isnan(ro(k))==1; continue; end;
                if bMinElevP(k)==0; continue; end;
                if isempty(bMinP{k})==1; continue; end;
                
                a1(k)=bMinElevP(k); % index for min cost within each Elevation level ~ similar to l index
                a2(k)=bMinP{k}(bMinElevP(k)); % index for min cost across ni Elevation levels ~ similar to m index
                
                if a2(k)==0; continue; end;
                if isnan(a2(k))==1; continue; end;
                
                aCOEPmin(k) = COEP{k}{a1(k)}(a2(k));
                aPPnetmin(k) = PPnet{k}{a1(k)}(a2(k));
                nPipeminP(k) = nPipeP{k}{a1(k)}(a2(k));
                
                
                if isempty(CostElementsP{k}{a1(k)})==1;
                    CostElementsPMin{k} = 0;
                    dfQPmin(k)= 0;
                    dfPLmin(k)= 0;
                    HeadPmin(k) = 0;
                    aPInletMin(k) = 0;
                    aPinlet_windowMin(k) = 0;
                    lat_Pin_min(k) = 0;
                    lon_Pin_min(k) = 0;
                    rPinMin(k) = 0;
                    cPinMin(k) = 0;
                    continue; end
                
                CostElementsPMin{k} = CostElementsP{k}{a1(k)}{a2(k)};
                dfQPmin(k)= QPinlet{k}{a1(k)}(a2(k));
                dfPLmin(k)= PL{k}{a1(k)}(a2(k));
                HeadPmin(k) = ZPinlet{k}{a1(k)}(a2(k)) - Zoutlets(k);
                HeadraceLPmin(k)= PL{k}{a1(k)}(a2(k)) - (ZPinlet{k}{a1(k)}(a2(k)) - Zoutlets(k));
                aPInletMin(k) = Pinlet{k}{a1(k)}(a2(k));
                aPinlet_windowMin(k) = Pinlet_win{k}{a1(k)}(a2(k));
                lat_Pin_min(k) = latP{k}{a1(k)}(a2(k));
                lon_Pin_min(k) = lonP{k}{a1(k)}(a2(k));
                rPinMin(k) = rem(aPInletMin(k)-1, nrw)+1;  %best inlet row
                cPinMin(k) = fix((aPInletMin(k)-1)/nrw)+1;  %best inlet col
                QDesignPinletMin(k) = QDesignPinlet{k}{a1(k)}(a2(k));
                QDesignLFPinletMin(k) = QDesignLFPinlet{k}{a1(k)}(a2(k));
                QDesignMeanPinletMin(k) = QDesignMeanPinlet{k}{a1(k)}(a2(k));
                ZPinletMin(k) = ZPinlet{k}{a1(k)}(a2(k));
                nPipePMin(k) = nPipeP{k}{a1(k)}(a2(k));
                OptInvPMin(k) = OptInvP{k}{a1(k)}(a2(k));
                PPMin(k) = PP{k}{a1(k)}(a2(k));
                OptSpecCapPMin(k) = OptSpecCapP{k}{a1(k)}(a2(k));
                accPmin(k) = accPinlet{k}{a1(k)}(a2(k));
                OpthfPMin(k) = OpthfP{k}{a1(k)}(a2(k));
                OptDPMin(k) = OptDP{k}{a1(k)}(a2(k));
                
            end
            
            %%%%Pcostmodel.m END
            %% Q-decreaser
            disp('Q decreaser for pipe systems');
            
            %%%%QPdecreaser.m START
            for k = 1:no
                fprintf('%d DP QP decreaser #%d of %d\n',nbasin,k,no)
                
                if isnan(ro(k))==1; continue; end;
                if isempty(aPinlet_windowMin(k)); continue; end
                if aPinlet_windowMin(k)==0; continue; end
                if isnan(aCOEPmin(k))==1; continue; end;
                
                clear r c
                [r,c] = ind2sub(size(Q_inlet_win{k}),aPinlet_windowMin(k));
                
                QQ = Qdecimater(r,c,inlet_win,Q_inlet_win{k},fdir_inlet_win{k},acc_inlet_win{k},aPinlet_windowMin(k));
                
                Q_inlet_win{k} = QQ;
                
            end
            
            %% Check
            % k=44;
            %
            % [rin,cin] = ind2sub(size(Q_inlet_win{k}),aPinlet_windowMin(k));
            %
            % figure(1);clf;imagesc(log(Q_inlet_win{k}));axis image; colormap(flipud(gray));
            % hold on
            % plot(inlet_win+1,inlet_win+1,'.r','markersize',20)
            % plot(cin,rin,'.b','markersize',20)
            % hold off
            
            %%%%QPdecreaser.m END
            %% Store info per Q decreaser loop
            disp('Store info per Q decreaser loop')
            
            PCOEend{ndl}=aCOEPmin; % $/kWh Dam-Pipe COE
            PPnetend{ndl}=aPPnetmin; % GWh Dam-Pipe energy potential
            Platminend{ndl} = lat_Pin_min;
            Plonminend{ndl}= lon_Pin_min;
            rPinMinend{ndl} = rPinMin;
            cPinMinend{ndl}= cPinMin;
            aPInletminEnd{ndl}=aPInletMin;
            aPinlet_windowMinEnd{ndl}=aPinlet_windowMin;
            dfQPminEnd{ndl}=dfQPmin;
            dfPLminEnd{ndl}=dfPLmin;
            CostElementsPMinEnd{ndl}=CostElementsPMin;
            Pinlet_winend{ndl}=Pinlet_win;
            HeadPminend{ndl}=HeadPmin;
            HeadraceLPminend{ndl}=HeadraceLPmin;
            QDesignPinletMinend{ndl}=QDesignPinletMin;
            QDesignLFPinletMinend{ndl}=QDesignLFPinletMin;
            QDesignMeanPinletMinend{ndl}=QDesignMeanPinletMin;
            dfPLminend{ndl}=dfPLmin;
            ZPinletMinend{ndl} = ZPinletMin;
            nPipePMinend{ndl} = nPipePMin;
            OptInvPMinend{ndl} = OptInvPMin;
            PPMinend{ndl} = PPMin;
            OptSpecCapPMinend{ndl} = OptSpecCapPMin;
            PSurfaceLake15sMinend{ndl} = PSurfaceLake15s;
            accPminend{ndl} = accPmin;
            OpthfPminend{ndl} = OpthfPMin;
            OptDminend{ndl} = OptDPMin; %Optimal pipe diameter
            
            %%
            aDPPnetEnd4  = horzcat(PPnetend{:})'; % GWh  Dam-Pipe systems Pnet
            fprintf('Total P potential: %0.0f GWh\n',sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
            
            %%%% PipeSystems.m END
            
        end
    end
    
    %% CLASSIC DAM SYSTEMS with RUN-OF-RIVER CHARACTERISTICS, LOADFACTOR < 100%
    if runRD==1
        disp('River dam scanner');
        
        %%%% riverdamscannerHR.m START
        %% Riverdam scanner to determine width and height in search radius
        
        counter=0;
        for  k = 1:numel(outlets)
            % k=32;
            fprintf('%d RD dam height-width #%d of %d\n',nbasin,k,numel(outlets))
            
            if isnan(ro(k))==1; continue; end;
            
            % First load high res dem
            
            RDwinsize_dem=30;
            %RDdemHR = dem_window_corners(latOut(k), lonOut(k), RDwinsize_dem, root_bil, continent_in); %Call in function to set up high res dem window
            RDdemHR = hs3s_window('con',latOut(k), lonOut(k), RDwinsize_dem, root_bil,continent_in);
            %RDdemHR = hs3s_window('dem',latOut(k), lonOut(k), RDwinsize_dem, root_bil,continent_in);
            RDdemHR = single(RDdemHR);
            
            
            % To control for return in demHR function
            isfile_outlets(k)=0;
            if RDdemHR==0;
                isfile_outlets(k)=1;
                ro(k) = NaN;
                co(k) = NaN;
                lonOut(k)=NaN;
                latOut(k)=NaN;
                dfhdamRD{k} = 0;
                dfldamRD{k} = 0;
                ZiCenterRD(k) = 0;
                continue;
            end;
            
            [nrRD,ncRD]=size(RDdemHR);
            
            % Second find lowest point in 15s highres window
            RDwin15s=2;
            
            firstrowRD_15s = max(1,RDwinsize_dem+1-RDwin15s);
            lastrowRD_15s  = min(nrRD,RDwinsize_dem+1+RDwin15s);
            firstcolRD_15s = max(1,RDwinsize_dem+1-RDwin15s);
            lastcolRD_15s  = min(ncRD,RDwinsize_dem+1+RDwin15s);
            
            RDdemHR_15s = RDdemHR(firstrowRD_15s:lastrowRD_15s,firstcolRD_15s:lastcolRD_15s); %Select small 15s window
            [RDr_mindemHR,RDc_mindemHR]=find(RDdemHR_15s==min(RDdemHR_15s(:))); %Find lowest point
            %Transformed focal point based on demHR_15s
            %RDrh_width = RDwinsize_dem+1; %r c of focus point of large high res dem
            %RDch_width = RDwinsize_dem+1;
            RDrh_width2 = RDwinsize_dem+1+(RDr_mindemHR(1)-RDwin15s-1); % New coordinates
            RDch_width2 = RDwinsize_dem+1+(RDc_mindemHR(1)-RDwin15s-1); % In case multiple lowest point, first option
            
            ZiCenterRD(k) = RDdemHR(RDrh_width2,RDch_width2); %Corrected Zfocus
            %figure(2);clf;imagesc(RDdemHR);axis image
            % %Check
            % figure(1);clf;imagesc(RDdemHR);axis image;colormap(jet);
            % figure(1);hold on;
            % plot(RDch_width,RDrh_width,'c.','markersize',20);
            % plot(RDch_width2,RDrh_width2,'r.','markersize',20);
            % hold off
            % figure(2);clf;imagesc(RDdemHR_15s);axis image;colormap(jet);
            % figure(2);hold on; plot(RDc_mindemHR,RDr_mindemHR,'r.','markersize',20); hold off
            
            %% Height difference in dem window used as dam height
            RDdemMax(k) = max(RDdemHR(:));
            
            if MiniHydro_special==1 & Qoutlets_design(k) < bighydro_cutoff
                MiniFlag(k)=1;
                if (RDdemMax(k)-ZiCenterRD(k)-1) > MiniDamMax;
                    RDZdiff(k) = MiniDamMax;
                    
                else
                    RDZdiff(k) = RDdemMax(k)-ZiCenterRD(k)-1;
                end
                
            else
                MiniFlag(k)=0;
                RDZdiff(k) = RDdemMax(k)-ZiCenterRD(k)-1;
            end
            
            dfhdamRD{k} = 1:RDZdiff(k)';
            
            
            %% Based on dam height we determine dam width
            
            if ASYMRD==0
                ldamRD=0;
                for wRD=1:RDZdiff(k)
                    % Determine Ztop based on hdam and Zfocus
                    clear RDZtop RDdx_w RDdy_w RDds_w RDhiground_w RDdhiground_w
                    RDZtop = ZiCenterRD(k) + wRD;
                    
                    % Build higround map
                    RDdx_w = repmat(abs(-RDwinsize_dem:RDwinsize_dem), [2*RDwinsize_dem+1, 1]);
                    RDdy_w = RDdx_w';
                    RDds_w = sqrt(RDdx_w.^2+RDdy_w.^2);
                    
                    % Determine minimum distance from Zfocus cell to higround
                    RDhiground_w = RDdemHR>RDZtop;
                    RDdhiground_w = RDds_w(RDhiground_w);
                    
                    ldamRD(wRD) = max(0.5,min(RDdhiground_w))*90*2;
                    
                    if isempty(ldamRD(wRD))==1;ldamRD(wRD)=0;end;
                    
                    %checkdfld
                    %fprintf('lenght = %f meters.\n',ldamRD(wRD))
                    %figure(1);clf;imagesc(RDhiground_w);axis image;colormap(jet);hold on
                    %plot(RDch_width2,RDrh_width2,'r.','markersize',20); hold off
                    
                end
                
                dfldamRD{k} = ldamRD;
                
            elseif ASYMRD==1
                
                [lrwidth{k} genwidth{k}]=AsymDamLenght(root_bil,continent_in,latOut(k),lonOut(k),RDwinsize_dem,RDZdiff(k),0);
                
                dfldamRD = lrwidth;
            end
            
        end
        
        
        %%%% riverdamscannerHR.m END
        
        
        %% River Dam population displacement
        disp('River dam population displacer');
        
        %%%% RDlake.m START
        %% RiverDam lake to calculate population displacement
        
        for k=1:numel(outlets)
            fprintf('%d RD Lake Pop-Landvalue #%d of %d\n',nbasin,k,numel(outlets))
            
            if isnan(ro(k))==1; continue; end;
            
            %Calculate upstream cells
            RDupstream = fastfindupstream_lim(acc,fdir,drow,dcol,outIdx(k));
            
            % define box around upstream cells to speed up
            ZupstreamRD = RDupstream==1;
            %ZupstreamRD = RDupstream==1 & Z > Z(ro(k),co(k)) & Z < (Z(ro(k),co(k))+dfhdamRD{k}(end));
            cs = find(sum(ZupstreamRD)>0); % columns
            rs = find(sum(ZupstreamRD')>0); % rows
            colwin=0;
            
            fcs = max(1,cs(1) - colwin);
            lcs = min(ncw,cs(end) + colwin);
            frs = max(1,rs(1) - colwin);
            lrs = min(nrw,rs(end) + colwin);
            Zupstream_win = Z(frs:lrs,fcs:lcs);
            RDupstream_win = RDupstream(frs:lrs,fcs:lcs);
            Popupstream_win = Pop(frs:lrs,fcs:lcs);
            LandValueupstream_win = LandValue(frs:lrs,fcs:lcs);
            
            %Determine which cells are flooded
            PopRDlake=0;
            MoneyRDlake=0;
            LandValueRDlakee=0;
            clear RDlake2 hd
            
            for j=1:numel(dfhdamRD{k})
                hd=single(dfhdamRD{k}(j));
                
                RDlake2 = RDupstream_win & Zupstream_win < (Zoutlets(k)-RDDepth(k)+OptDH(k)); % Based on 15s DEM map, minus river depth + dam height based on 3s DEM
                
                %Surface
                RDlakeSurface{k}(j) = sum(RDlake2(:)==1)* (450*(450*cosd(latOut(k)))); % m2
                
                %Volume
                Z_upstream = Zupstream_win(RDupstream_win);
                dz_h = Zoutlets(k)-RDDepth(k)+hd-Z_upstream;
                dz_h = max(0, dz_h);
                RDVolumeLake{k}(j)  = sum(dz_h * 450 * (450*cosd(latOut(k)))); % m3
                
                %Population
                RDlakeidx = find(RDlake2);
                PopRDlake(j) = sum(Popupstream_win(RDlakeidx));
                
                %Land value check with runsettings
                if landval_set==1
                    LandValueRDlakee(j) = sum(LandValueupstream_win(RDlakeidx))*0.45^2 * 9; %LandValue is $/km2. Area of 15s cell is 450m^2 x 9.09 (infinite discounted yearly revenue with discount factor 10%, see perpetuity_check.m)
                else
                    LandValueRDlakee(j) = 0;
                end
                
                %% Find GDP values
                %GDPSwiss = 55000; % $/cap/year
                %GDPUSA = 54000; % $/cap/year
                [GDPr,~] = find(ISOGDP(:,1)==RDCountry_id(k));
                GDPpc = ISOGDP(GDPr,3);
                if isempty(GDPr)==1; GDPpc = mean(ISOGDP(:,3)); end %if ISO value is not found, world average
                
                MoneyRDlake(j) = 2 * GDPpc * PopRDlake(j); %
            end
            
            PopDisplaced{k}=PopRDlake;
            PopCost{k}=MoneyRDlake;
            %RDLake_map{k}=RDlake2;
            LandValueRDlake{k}=LandValueRDlakee;
            
        end
        
        %% Check lake
        
        % figure(1);clf;
        % img1= truecolorsc(Zupstream_win,flipud(gray));
        % % img2 = burnmask(img1, ~RDLake_map{14}{40});
        % img2 = burnmask(img1, ~RDlake2);
        % %             img2 = burnmask(img1, ~RDLake_map{k}{200});
        % image(img2); axis image;
        % title('Discharge','FontSize',14);
        % %     axis([75 125 75 125]);
        % set(gca,'xtick',[],'ytick',[])
        % hold on;
        % % plot(co(23),ro(23),'r.','markersize',15);
        % hold off
        
        %%
        % [nr nc] = size(Zupstream_win)
        % for r=1:nr
        %     for c=1:nc
        %         if RDupstream_win(r,c)==1
        %             Zup(r,c) = Zupstream_win(r,c);
        %         else
        %             Zup(r,c) = NaN;
        %         end
        %     end
        % end
        % figure(3);clf;imagesc(Zup); axis image
        %
        % %%
        % for r=1:nr
        %     for c=1:nc
        %         Zupm(r,c) = 710 - Zup(r,c);
        %     end
        % end
        % figure(4);clf;imagesc(Zupm); axis image
        %
        % %%
        % Zupmm = max(0, Zupm);
        % figure(5);clf;imagesc(Zupmm); axis image
        
        %%%% RDlake.m END
        
        %% River Dam scanner
        disp('River dam Cost model');
        
        %%%%CostModel_RiverDam.m START
        %% Costmodel_RiverDam
        % Calculate cost and potential of river dam
        
        COETotRD = 0;
        RDPnet=0;
        clear DH
        for  k = 1:no
            
            if isnan(ro(k))==1;
                COETotRD(k)=0;
                RDPnet(k)=0;
                OptDH(k)=0;
                OptDL(k)=0;
                RDCostElements{k}=0;
                continue;
            end; % Skip empty arrays
            
            if numel(dfldamRD{k})==0;
                COETotRD(k)=0;
                RDPnet(k)=0;
                OptDH(k)=0;
                OptDL(k)=0;
                RDCostElements{k}=0;
                continue;
            end; % Skip empty arrays
            
            if numel(dfhdamRD{k})==0;
                COETotRD(k)=0;
                RDPnet(k)=0;
                OptDH(k)=0;
                OptDL(k)=0;
                RDCostElements{k}=0;
                continue;
            end; % Skip empty arrays
            
            DH= single(dfhdamRD{k});
            
            if MiniHydro_special==1 & Qoutlets_design(k) < bighydro_cutoff
                [COETotRD(k) RDPnet(k) OptP(k) OptDH(k) OptDL(k) RDCostElements{k}]=MiniCostmodel(dfldamRD{k},DH,PopCost{k},Qoutlets_design(k),Qoutlets_design_LF(k),DisOutlet(k),LandValueRDlake{k},RDDepth(k),nbasin,k,cost_constr,cost_lim);
            else
                [COETotRD(k) RDPnet(k) OptP(k) OptDH(k) OptDL(k) RDCostElements{k} OptInv(k) RDP(k) OptPop(k) OptLV(k) OptSpecCap(k)]=RDcostmodel(dfldamRD{k},DH,PopCost{k},Qoutlets_design(k),Qoutlets_design_LF(k),DisOutlet(k),LandValueRDlake{k},RDDepth(k),nbasin,k,cost_constr,cost_lim,Quakerate(k));
            end
            
            PopDisplacedOpt(k) = PopDisplaced{k}(OptDH(k));
        end
        
        
        %% To remove zeros
        for k =1:no
            if COETotRD(k)==0; COETotRD(k)=NaN; end;
            if RDPnet(k)==0; RDPnet(k)=NaN;end;
        end
        
        COETotRD = COETotRD';
        RDPnet = RDPnet'; % kWh
        RDPnet = RDPnet*1e-6; % GWh
        
        
        
        %%%%CostModel_RiverDam.m END
        
        %% River dam selector
        disp('Dam selector');
        if runMini==0
            %%%% Dam_selector.m START
            %% Dam selector
            % Routine to de-selects dams according to a certain priority
            
            if sum(isnan(COETotRD))==numel(ro)
                RDlake_Opt=0;
                CheapestDam=0;
                DeselectedSites=0;
                DeselectedSites_unSort=0;
                COETotRDs=0;
                RDPnets=0;
            end
            
            %% Collecting all dam locations DP, P and RD
            
            if runDP==1
                %Zeros are NaN
                aInletminEnd{1}(aInletminEnd{1}(:)==0)=NaN;
                aInletminEnd{2}(aInletminEnd{2}(:)==0)=NaN;
                outlets(outlets(:)==0)=NaN;
                
                DPlocs = horzcat(aInletminEnd{:})';
                Plocs  = horzcat(aPInletminEnd{:})';
                RDlocs = outlets';
                
                Damlocs = [DPlocs;Plocs;RDlocs]; %In terms of basin indices
                
                % Vector to keep track of type of systems
                SysID  = zeros((numel(RDlocs)+numel(DPlocs)+numel(Plocs)),1);
                SysID(1:numel(DPlocs))=1;                                      %DP systems
                SysID((numel(DPlocs)+1):(numel(DPlocs)+numel(Plocs)))=2;       %P systems
                SysID((numel(DPlocs)+numel(Plocs)+1):end)=3;                   %RD systems
                
                % Rows and columns
                ros = [horzcat(rinMinend{:})'; horzcat(rPinMinend{:})'; ro];
                cos = [horzcat(cinMinend{:})'; horzcat(cPinMinend{:})'; co];
                
                lats = [horzcat(latminend{:})'; horzcat(Platminend{:})'; latOut];
                lons = [horzcat(latminend{:})'; horzcat(Plonminend{:})'; lonOut];
                
                % %%
                % figure(1);clf;imagesc(log(Q));colormap(flipud(gray));axis image;
                % hold on
                % plot(co(find(~isnan(COETotRD))),ro(find(~isnan(COETotRD))),'.r','markersize',20)
                % plot(cinMinend{1}(find(~isnan(COEend{1}))),rinMinend{1}(find(~isnan(COEend{1}))),'.b','markersize',20)
                % plot(cPinMinend{1}(find(~isnan(PCOEend{1}))),rPinMinend{1}(find(~isnan(PCOEend{1}))),'.g','markersize',20)
                % hold off
                
            else
                %Zeros are NaN
                outlets(outlets(:)==0)=NaN;
                
                Plocs  = horzcat(aPInletminEnd{:})';
                RDlocs = outlets';
                
                Damlocs = [Plocs;RDlocs]; %In terms of basin indices
                Grandlocs = GrandIdx; %Existing locations
                
                % Vector to keep track of type of systems
                SysID  = zeros((numel(RDlocs)+numel(Plocs)),1);
                SysID((1:numel(Plocs)))=1;       %P systems
                SysID((numel(Plocs)+1):end)=2;   %RD systems
                
                % Rows and columns
                ros = [horzcat(rPinMinend{:})'; ro];
                cos = [horzcat(cPinMinend{:})'; co];
                
                lats = [horzcat(Platminend{:})'; latOut'];
                lons = [horzcat(Plonminend{:})'; lonOut'];
                
            end
            
            %% Show
            % [rRD,cRD]=ind2sub(size(acc),RDlocs);
            % % [rDP,cDP]=ind2sub(size(acc),DPlocs);
            % [rP,cP]=ind2sub(size(acc),Plocs);
            %
            % figure(1);clf;imagesc(log(Q));axis image;colormap(flipud(gray));
            % hold on
            % % plot(cRD,rRD,'r.','markersize',10)
            % % plot(cDP,rDP,'b.','markersize',20)
            % % plot(cP,rP,'g.','markersize',10)
            % plot(coss,ross,'g.','markersize',20)
            % plot(c_damst,r_damst,'b.','markersize',20)
            % hold off
            %% Collecting COEs and Pnets
            
            if runDP==1
                % COEs
                COEendDP = horzcat(COEend{:});
                COEendDP(COEendDP(:)==0)=NaN;
                
                COEendP = horzcat(PCOEend{:});
                COEendP(COEendP(:)==0)=NaN;
                
                COEAll = [COEendDP';COEendP';COETotRD];
                
                % Pnets
                PnetendDP = horzcat(DPPnetend{:});
                PnetendDP(PnetendDP(:)==0)=NaN;
                
                PnetendP = horzcat(PPnetend{:});
                PnetendP(PnetendP(:)==0)=NaN;
                
                PnetAll = [PnetendDP';PnetendP';RDPnet];
                
            else
                % COEs
                COEendP = horzcat(PCOEend{:});
                COEendP(COEendP(:)==0)=NaN;
                
                COEAll = [COEendP';COETotRD];
                
                % Pnets
                PnetendP = horzcat(PPnetend{:});
                PnetendP(PnetendP(:)==0)=NaN;
                
                PnetAll = [PnetendP';RDPnet];
                
            end
            
            %% First, recalculate lake based on optimal dam height and check which lake floods which dams
            %% DP-systems
            
            if runDP==1
                if sum(isnan(COEendDP))~=numel(COEendDP)
                    for ndd=1:2
                        for k=1:no
                            fprintf('Recalculating DP lakes outlet #%d of %d\n',k,no);
                            
                            clear DPupstream a OutletDeselect
                            if isnan(ro(k))==1;
                                DeselectedSites_unSortDP{ndd}{k}=0;
                                continue; end;
                            if isnan(aInletminEnd{ndd}(k))==1;
                                DeselectedSites_unSortDP{ndd}{k}=0;
                                continue; end;
                            
                            DPupstream = fastfindupstream_lim(acc,fdir,drow,dcol,aInletminEnd{ndd}(k));
                            
                            DPlake_Opt = DPupstream & Z < (Z(aInletminEnd{ndd}(k))-DPDepth{ndd}(k)+ahdammin{ndd}(k));
                            
                            LakeIdx = find(DPlake_Opt);
                            a = ismember(Damlocs,LakeIdx);
                            OutletDeselect = find(a);
                            kc=k;
                            if ndd==2; kc=k+no;end % Making sure it doesn't select itself in the second round
                            DeselectedSites_unSortDP{ndd}{k} = OutletDeselect(OutletDeselect~=kc); % Which dams are flooded (unsorted)
                            
                            clear a
                            a = ismember(Grandlocs,LakeIdx);
                            GrandOutletDeselectDO = find(a);
                            
                        end
                    end
                end
            end
            
            %% Show stuation before deselection
            
            % for i=1:numel(DPlake_Opt{1})
            %     a(i)=isempty(DPlake_Opt{1}{i});
            % end
            % b=find(~a);
            %
            % %
            % figure(1);clf;
            % img{1}= truecolorsc(DPlake_Opt{1}{b(1)},flipud(gray));
            % %
            % for i=1:numel(b)
            %     if i==1; continue; end
            %         img{i} = burnmask(img{i-1}, ~DPlake_Opt{1}{b(i)});
            % end
            % image(img{end}); axis image;
            % %
            % hold on
            % plot(cinMinend{1}(:),rinMinend{1}(:),'.r','markersize',15)
            % plot(cPinMinend{1}(:),rPinMinend{1}(:),'.g','markersize',15)
            % plot(cinMinend{2}(:),rinMinend{2}(:),'.r','markersize',15)
            % plot(cPinMinend{2}(:),rPinMinend{2}(:),'.g','markersize',15)
            % plot(co,ro,'.b','markersize',15)
            % hold off
            
            
            %% RD-systems
            fprintf('# of valid RD-systems: %d of %d\n',sum(isnan(COETotRD)),no);
            if sum(isnan(COETotRD))~=numel(ro)
                for k=1:no
                    fprintf('Recalculating lakes outlet #%d of %d\n',k,no);
                    
                    clear a OutletDeselect Z_upstream RDupstream dz_h RDlake_Opt
                    if isnan(ro(k))==1;
                        DeselectedSites_unSortRD{k}=0;
                        DeselectedGrandSites_unSortRD{k}=0;
                        continue; end;
                    if isnan(COETotRD(k))==1;
                        DeselectedSites_unSortRD{k}=0;
                        DeselectedGrandSites_unSortRD{k}=0;
                        continue; end;
                    
                    RDupstream = fastfindupstream_lim(acc,fdir,drow,dcol,outlets(k));
                    
                    RDlake_Opt = RDupstream & Z < (Zoutlets(k)-RDDepth(k)+OptDH(k));
                    
                    %For potential dams
                    LakeIdx = find(RDlake_Opt);
                    a = ismember(Damlocs,LakeIdx);
                    OutletDeselect = find(a);
                    if runDP==1
                        DeselectedSites_unSortRD{k} = OutletDeselect(OutletDeselect~=(numel(DPlocs)+numel(Plocs)+k)); % Which dams are flooded (unsorted)
                    else
                        DeselectedSites_unSortRD{k} = OutletDeselect(OutletDeselect~=(numel(Plocs)+k)); % Which dams are flooded (unsorted)
                    end
                    
                    %For existing dams
                    clear a
                    a = ismember(Grandlocs,LakeIdx);
                    GrandOutletDeselectRD = find(a); %
                    DeselectedGrandSites_unSortRD{k} = GrandOutletDeselectRD; % Which Grand dams are flooded (unsorted)
                    
                    %% Lake volume calculator m3
                    Z_upstream = Z(RDupstream);
                    dz_h = Zoutlets(k)-RDDepth(k)+OptDH(k)-Z_upstream;
                    dz_h = max(0, dz_h);
                    RDVolumeLake15s(k)  = sum(dz_h * 450 * (450*cosd(latOut(k)))); %Volume lake m3
                    RDSurfaceLake15s(k) = max(1,sum(RDlake_Opt(:))) * 450 * (450*cosd(latOut(k))); %Surface Reservoir m2
                end
            end
            
            %% Collecting all deselected DamLocs
            %clear DeselectedSites CheapestDam DeselectedSites_final
            
            if runDP==1
                DeselectedSites_unSortDPP = horzcat(DeselectedSites_unSortDP{:});
                
                for k = 1:numel(DeselectedSites_unSortDPP); DeselectedSites_unSortP{k}=[]; end %Creating empty P array (floods nothing)
                
                DeselectedSites_unSort = [DeselectedSites_unSortDPP,DeselectedSites_unSortP,DeselectedSites_unSortRD];
                
            else
                for k = 1:(nd*numel(DeselectedSites_unSortRD)); DeselectedSites_unSortP{k}=[]; end %Creating empty P array (floods nothing)
                
                DeselectedSites_unSort = [DeselectedSites_unSortP,DeselectedSites_unSortRD];
                DeselectedGrandSites_unSort = [DeselectedSites_unSortP,DeselectedGrandSites_unSortRD];
                
                % Lake surfaces
                LakeSurfacesAll = [horzcat(PSurfaceLake15sMinend{:})'; RDSurfaceLake15s'];
                
                % Flow accumulations
                accAll = [horzcat(accPminend{:})'; Accoutlets'];
                
                % Basin IDs
                BIDall = zeros(size(accAll));
                BIDall(:) = nbasin;
                
                % Continent IDs
                CIDall = zeros(size(accAll));
                CIDall(:) = CID;
            end
            
            %% Show stuation before deselection
            % figure(1);clf;
            %
            % img{1}= truecolorsc(RDlake_Opt{1},flipud(gray));
            %
            % for i=1:numel(RDlake_Opt)
            %     if i==1; continue; end
            %
            %     img{i} = burnmask(img{i-1}, ~RDlake_Opt{i});
            % end
            % image(img{end}); axis image;
            %
            % hold on
            % plot(co(1:i),ro(1:i),'.r','markersize',15)
            % hold off
            
            %% Second, sorting through indexing
            
            if sum(isnan(COETotRD))~=numel(ro)
                if exist('DeselectedSites_unSort')
                    
                    outletsRel=find(~isnan(COEAll)); %Find relevant outlet idx
                    
                    %Cost Priority
                    minCOE = min(COEAll);
                    COEIndexing = 1./(COEAll/minCOE); %Lowest COE gets 1
                    if sum(isnan(COEIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in COE deselector');return;end
                    if sum(isinf(COEIndexing(outletsRel))) > 0; returnID=1;disp('Infs in COE deselector');return;end
                    
                    %Power priority
                    maxPnet = max(PnetAll);
                    PnetIndexing = PnetAll/maxPnet; %Highest Pnet gets 1
                    if sum(isnan(PnetIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in Pnet deselector');return;end
                    if sum(isinf(PnetIndexing(outletsRel))) > 0; returnID=1;disp('Infs in Pnet deselector');return;end
                    
                    %Cost per power priority
                    COEperPnet = COEAll./PnetAll;
                    minCOEperPnet = min(COEperPnet);
                    COEperPnetIndexing = 1./(COEperPnet/minCOEperPnet); %Lowest COEperPnet gets 1
                    if sum(isnan(COEperPnetIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in COEPnet deselector');return;end
                    if sum(isinf(COEperPnetIndexing(outletsRel))) > 0; returnID=1;disp('Infs in COEPnet deselector');return;end
                    
                    %Lake surface priority
                    MaxSurfaceLake = max(LakeSurfacesAll);
                    LakeIndexing = 1./(LakeSurfacesAll/MaxSurfaceLake); %Lowest surface gets 1
                    if sum(isnan(LakeIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in Lake deselector');return;end
                    if sum(isinf(LakeIndexing(outletsRel))) > 0; returnID=1;disp('Infs in Lake deselector');return;end
                    
                    %Flow accumulation priority
                    MaxaccAll = max(accAll);
                    AccIndexing = 1./(accAll/MaxaccAll); %Lowest flow acc gets 1
                    if sum(isnan(AccIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in acc deselector');return;end
                    if sum(isinf(AccIndexing(outletsRel))) > 0; returnID=1;disp('Infs in acc deselector');return;end
                    
                    WeightedCOEPnet = (COEIndexing*betaC)+(PnetIndexing*betaP)+(COEperPnetIndexing*betaCP)+(LakeIndexing*betaL)+(AccIndexing*betaA);
                    
                    %Sorting
                    [~, sortIdx] = sort(WeightedCOEPnet(outletsRel),'descend');
                    
                    CheapestDam = outletsRel(sortIdx);
                    DeselectedSites = DeselectedSites_unSort(outletsRel(sortIdx)); %Sorted flooded dams
                    
                else
                    RDlake_Opt=0;
                    CheapestDam=0;
                    DeselectedSites=0;
                    COEAll=0;
                    PnetAll=0;
                end
            end
            
            %% Third, deselection starting with cheapest
            if iscell(DeselectedSites)
                CheapestDam2=CheapestDam;
                for i=1:numel(CheapestDam2)
                    if sum(ismember(CheapestDam2(i+1:end), DeselectedSites{i}))>=1;
                        idxdeselect{i} = find(ismember(CheapestDam2(i+1:end), DeselectedSites{i})) + i;
                        CheapestDam2(idxdeselect{i})=1;
                        continue; end
                end
            end
            
            %% Fourth, deselection starting with most expensive
            if iscell(DeselectedSites)
                CheapestDam3=CheapestDam2;
                for i=0:numel(CheapestDam3)-1
                    if i==numel(CheapestDam3)-1; continue; end
                    if sum(ismember(CheapestDam3(:), DeselectedSites{numel(CheapestDam3)-i})) >=1
                        CheapestDam3(numel(CheapestDam3)-i)=1;
                    end
                end
            end
            
            %% Final check, none of the dams should be flooded
            if iscell(DeselectedSites)
                for i=1:numel(CheapestDam)
                    if CheapestDam3(i)==1; continue; end
                    DeselectedSites_final{i} = DeselectedSites{i};
                end
                
                for i=1:numel(DeselectedSites_final)
                    if isrow(DeselectedSites_final{i})==0; DeselectedSites_final{i}=DeselectedSites_final{i}'; end
                end
                
                
                if exist('DeselectedSites_final')
                    DeselectedSites_Final_vector = horzcat(DeselectedSites_final{:})';
                    FloodCheck = sum(ismember(CheapestDam3(:),DeselectedSites_Final_vector(:)));
                    fprintf('Flood check = %0.0f\n',FloodCheck);
                end
            end
            
            %% Check which index gives problem
            
            % for i=1:numel(CheapestDam3)
            %     x=sum(find(DeselectedSites_final{i}==899));
            %     if x>=1
            %         i
            %     end;
            % end
            
            %% New COETotRD and RDPnet
            COEAlls = zeros(size(COEAll));
            COEAlls(COEAlls==0)=NaN;
            PnetAlls = zeros(size(PnetAll));
            PnetAlls(PnetAlls==0)=NaN;
            ross = zeros(size(ros));
            ross(ross==0)=NaN;
            coss = zeros(size(cos));
            coss(coss==0)=NaN;
            SysIDs = zeros(size(SysID));
            SysIDs(SysIDs==0)=NaN;
            latss = zeros(size(lats));
            latss(latss==0)=NaN;
            lonss = zeros(size(lons));
            lonss(lonss==0)=NaN;
            LakeSurfacesAlls = zeros(size(LakeSurfacesAll));
            LakeSurfacesAlls(LakeSurfacesAlls==0)=NaN;
            accAlls = zeros(size(accAll));
            accAlls(accAlls==0)=NaN;
            BIDs = zeros(size(BIDall));
            BIDs(BIDs==0)=nbasin;
            CIDs = zeros(size(CIDall));
            CIDs(CIDs==0)=CID;
            
            if iscell(DeselectedSites)
                for i=1:numel(CheapestDam3)
                    if CheapestDam3(i)>1
                        COEAlls(CheapestDam3(i))=COEAll(CheapestDam3(i));
                        PnetAlls(CheapestDam3(i))=PnetAll(CheapestDam3(i));
                        
                        SysIDs(CheapestDam3(i))=SysID(CheapestDam3(i));
                        ross(CheapestDam3(i))=ros(CheapestDam3(i));
                        coss(CheapestDam3(i))=cos(CheapestDam3(i));
                        
                        latss(CheapestDam3(i))=lats(CheapestDam3(i));
                        lonss(CheapestDam3(i))=lons(CheapestDam3(i));
                        
                        LakeSurfacesAlls(CheapestDam3(i))=LakeSurfacesAll(CheapestDam3(i));
                        
                        accAlls(CheapestDam3(i))=accAll(CheapestDam3(i));
                        BIDs(CheapestDam3(i))=BIDall(CheapestDam3(i));
                        CIDs(CheapestDam3(i))=CIDall(CheapestDam3(i));
                    end
                end
                
                fprintf('Total potential before deselection: %0.0f GWh\n',sum(PnetAll(~isnan(PnetAll))));
                fprintf('Total potential after deselection: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
                
                %Deselecting dams that flood existing Grand dams
                for i=1:numel(DeselectedGrandSites_unSort)
                    if isempty(DeselectedGrandSites_unSort{i})==1; continue; end
                    for j=1:numel(DeselectedGrandSites_unSort{i})
                        %fprintf('%d %d. %0.2f\n',i,DeselectedGrandSites_unSort{i}(j),COEAlls(i))
                        COEAlls(i)=NaN;
                        PnetAlls(i)=NaN;
                        SysIDs(i)=NaN;
                        ross(i)=NaN;
                        coss(i)=NaN;
                        latss(i)=NaN;
                        lonss(i)=NaN;
                        LakeSurfacesAlls(i)=NaN;
                        accAlls(i)=NaN;
                        BIDs(i)=NaN;
                        CIDs(i)=NaN;
                    end
                end
                fprintf('Total potential after deselecting Grand: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
                
                DPIDs = numel(find(SysIDs==1));
                PIDs = numel(find(SysIDs==2));
                RDIDs = numel(find(SysIDs==3));
            end
            
            %% Show dams and lake
            % figure(2)
            % clear img
            % img{1}= truecolorsc(RDlake_Opt{CheapestDam3(1)},flipud(gray));
            % for i=1:numel(CheapestDam3)
            %     if i==1; continue; end
            %     if CheapestDam3(i)==1;
            %         img{i} = img{i-1};
            %         continue;
            %     end
            %
            %     img{i} = burnmask(img{i-1}, ~RDlake_Opt{CheapestDam3(i)});
            % end
            % %
            % image(img{end}); axis image;
            % title('Second selection')
            % hold on
            % plot(co(CheapestDam3),ro(CheapestDam3),'.r','markersize',15)
            % hold off
            %
            % %% Show top-tot-teen graph
            % % Dams can still overlap when they're in different streams
            % indices_select = CheapestDam(find(~isnan(ro(CheapestDam3))));
            % figure(5); clf;
            % hold on;
            % for i=1:numel(indices_select)
            %     zbot = Z(ro(indices_select(i)),co(indices_select(i)));
            %     ztop = zbot + OptDH(indices_select(i));
            %     carea = acc(ro(indices_select(i)),co(indices_select(i)));
            %     x = [log(carea) log(carea)];
            %     y = [zbot ztop];
            %     plot(x,y, 'linewidth',2,'color',rand(1,3));
            % end;
            % hold off;
            % xlabel('log(area)');
            % ylabel('Elevation');
            
            % [lo,la] = setltln(acc, Rw, r_damst(10), c_damst(10)); %Coordinates outlets
            
            %%%% Dam_selector.m END
        end
        
        if returnID; continue; end
        
    end
    
    %% Save and sum
    disp('Saving energy data')
    
    %     matpath = fullfile(root, sprintf('output\\%s', continent_in));
    %     if ~isdir(matpath)
    %         mkdir(matpath);
    %     end
    
    %% Save Energy data
    
    if runDP==1 & runP==1 & runRD==1 & runMini==0
        save(matfile,'-v7.3' ...
            ,'COEend','DPPnetend','lat','lon','DPRegion_id','DPCountry_id','rin','cin','rinMinend','cinMinend','lonminend','latminend','aInletminEnd','ainlet_windowMinEnd','ahdammin','aldammin','CostElementsMinend','ainlet_windowMin','rinMin','cinMin','DPDepth' ... % DP vars
            ,'COETotRD','RDPnet','Qoutlets_design','Qoutlets_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets','DisOutlet','OptSpecCap' ... % RD vars
            ,'aPinlet_windowMinEnd','PCOEend','PPnetend','Platminend','Plonminend','rPinMinend','cPinMinend','aPInletminEnd','CostElementsPMinEnd','aPinlet_windowMin','OptInvP','OptInvPMinend','PP','PPMinend','Pinlet_winend','dfQPminEnd','dfPLminEnd','HeadPminend','QDesignPinletMinend','QDesignLFPinletMinend','QDesignMeanPinletMinend','ZPinletMinend','nPipePMinend' ... % P vars
            ,'COEAll','PnetAll','SysID','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss' ... % General vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','inlet_window','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    elseif runDP==0 & runP==1 & runRD==1 & runMini==0 %Standard output
        save(matfile,'-v7.3' ...
            ,'COETotRD','RDPnet','Qoutlets_design','Qoutlets_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','OptPop','OptLV','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets','DisOutlet','OptSpecCap','PopDisplacedOpt','RDlakeSurface','RDVolumeLake' ... % RD vars
            ,'aPinlet_windowMinEnd','PCOEend','PPnetend','Platminend','Plonminend','rPinMinend','cPinMinend','aPInletminEnd','CostElementsPMinEnd','aPinlet_windowMin','OptInvP','OptInvPMinend','PP','PPMinend','Pinlet_winend','dfQPminEnd','dfPLminEnd','HeadPminend','HeadraceLPminend','QDesignPinletMinend','QDesignLFPinletMinend','QDesignMeanPinletMinend','ZPinletMinend','nPipePMinend','OptSpecCapPMinend','rPin','cPin','OpthfPminend', 'OptDminend'... % P vars
            ,'COEAll','PnetAll','SysID','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss','LakeSurfacesAlls','accAlls' ... % General vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
        save(matfileCOEPOT,'-v7.3' ...
            ,'COEAlls','PnetAlls','SysIDs','RDDepth','RDCountry_id','RDRegion_id','latss','lonss'); %Just COE and Pnet
    elseif runDP==1 & runP==0 & runRD==1 && runMini==0
        save(matfile,'-v7.3' ...
            ,'COEend','DPPnetend','lat','lon','DPRegion_id','DPCountry_id','rin','cin','rinMinend','cinMinend','lonminend','latminend','aInletminEnd','ainlet_windowMinEnd','ahdammin','aldammin','CostElementsMinend','ainlet_windowMin','rinMin','cinMin','DPDepth' ... % DP vars
            ,'COETotRD','RDPnet','Qoutlets_design','Qoutlets_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets','DisOutlet','OptSpecCap' ... % RD vars
            ,'COEAll','PnetAll','SysID','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss','BIDs','CIDs' ... % General vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','inlet_window','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    elseif runDP==0 & runP==0 & runRD==1 & runMini==0
        save(matfile,'-v7.3' ...
            ,'COETotRD','RDPnet','Qoutlets_design','Qoutlets_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets','DisOutlet','OptSpecCap' ... % RD vars
            ,'COEAll','PnetAll','SysID','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss' ... % General vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','inlet_window','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    elseif runDP==0 & runP==1 & runRD==0 & runMini==0
        save(matfile,'-v7.3' ...
            ,'latOut','lonOut','ro','co','aPinlet_windowMinEnd','PCOEend','PPnetend','Platminend','Plonminend','rPinMinend','cPinMinend','aPInletminEnd','CostElementsPMinEnd','aPinlet_windowMin','OptInvP','PP','Pinlet_winend','dfQPminEnd','dfPLminEnd','HeadPminend','QDesignPinletMinend','QDesignLFPinletMinend' ... % P vars
            ,'nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    elseif runDP==0 & runP==0 & runRD==1 & runMini==1
        save(matfile,'-v7.3' ...
            ,'COETotRD','RDPnet','Qoutlets_design','Qoutlets_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','OptSpecCap' ... % RD vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    end
    
    %% Report
    if runMini==0 && runP==0
        fprintf('Tot potential: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
    elseif runP==1
        fprintf('Tot PnetAlls potential: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
        fprintf('Tot P potential: %0.0f GWh\n',sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
        fprintf('Tot RD potential: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet))));
        fprintf('Tot RD + P potential: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet)))+sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
    else
        fprintf('Tot RD potential: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet))));
    end
    
    %     catch ME %Error handling
    %
    %        fprintf('\nError in basin %d:\n\n%s\n\n',nbasin, ME.message);
    %
    %        errormessage = sprintf('%s',ME.message);
    %
    %        errpath = fullfile(root, sprintf('output\\%s\\%s',version, continent_in));
    %        if ~isdir(errpath)
    %            mkdir(errpath);
    %        end
    %
    %        txtfile = fullfile(root, sprintf('output\\%s\\%s\\Error_%d.txt',version, continent_in, nbasin));
    %        fileID = fopen(txtfile,'w');
    %        fprintf(fileID,'%s',errormessage);
    %        fclose(fileID);
    
    %%   end
end
toc
diary off
if saveProfile
    profsave(profile('info') ,fullfile(matpath,'profile_AllInOne'))
end

fprintf('\n HYDRUS01 END TEST RUN for do=%d:  %s\n',do,datetime('now'));
