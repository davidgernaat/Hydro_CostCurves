function NoMainstrm = FirstMainstrmDam(Q,adir,acc,r_damst,c_damst,GrandIdx)

%% Find mainstream
[~,MainIdx] = max(acc(:));
mainstrm = find_mainstream(Q,MainIdx,adir);

b=find(ismember(GrandIdx,mainstrm));

if isempty(b)==1; NoMainstrm=1; return; end

%% Find idx closest to outlet
for i=1:numel(b)
    GrandAcc(i) = acc(r_damst(b(i)),c_damst(b(i)));
end

[~,MaxAccIdx]=max(GrandAcc);

%% Deselect all mainstream cells prior to first dam
MainstrmAccs = acc(mainstrm);
MainstrmAccs(MainstrmAccs>GrandAcc(MaxAccIdx))=NaN;

MainstrmDes = find(isnan(MainstrmAccs(:))); %Deselect these cells from mainstream

NoMainstrm = mainstrm(MainstrmDes);

%%
% figure(1);clf;
% a=zeros(size(Q),'single');
% a(fdir(:)==9)=-1;
% a(NoMainstrm)=1;
% imagesc(a);colormap(gray);axis image;hold on
% plot(c_damst,r_damst,'.r','markersize',20)
% plot(c_damst(b),r_damst(b),'.b','markersize',20)
% plot(c_damst(b(MinDisIdx)),r_damst(b(MinDisIdx)),'.c','markersize',20)
% hold off

end