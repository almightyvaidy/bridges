clear;
clc;
close all;
tic
thisrun=string(datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm-ss'));
outfolder=strcat("C:\Users\vaidy\Google Drive\bridges\", thisrun,"\");
mkdir(outfolder)
cityvalue=4;
stateval=0;
switch stateval
    case 0 
        xd0=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_0_final - Copy (2).csv');toc
        xd1=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_1_final - Copy (2).csv');toc
        xd2=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_2_final - Copy (2).csv');toc
        xd3=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_3_final - Copy (2).csv');toc
        xd4=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_4_final - Copy (2).csv');toc
        xd5=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_5_final - Copy (2).csv');toc
        xd6=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_6_final - Copy (2).csv');toc
        xd7=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_7_final - Copy (2).csv');toc
        xd8=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_8_final - Copy (2).csv');toc
        xd9=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_9_final - Copy (2).csv');toc
        xd10=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_10_final - Copy (2).csv');toc
        xd11=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_11_final - Copy (2).csv');toc
        xd12=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_12_final - Copy (2).csv');toc
        xd13=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_13_final - Copy (2).csv');toc
        xd14=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_14_final - Copy (2).csv');toc
        xd15=readtable('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\Texas\texas_15_final - Copy (2).csv');toc
        xdbig=[xd0;xd1;xd2;xd3;xd4;xd5;xd6;xd7;xd8;xd9;xd10;xd11;xd12;xd13;xd14;xd15];toc;
    case 1
end
%xdbig(xdbig.bridge== 'F') = 0;
%xdbig{xdbig.bridge == 'T'} = 1;
isbridge=strcmp(xdbig.bridge,'T');
istunnel=strcmp(xdbig.tunnel,'T');
xbig=[table2array(xdbig(:,1)) table2array(xdbig(:,2)) isbridge istunnel table2array(xdbig(:,5)) table2array(xdbig(:,6)) table2array(xdbig(:,7)) table2array(xdbig(:,8)) ];toc;
 city=[ -95.3103 29.7752;-97.7341 30.2849; -96.6745	32.7100;  -98.525262	29.59002]
% 
%  lonsandlats=[-97.7394	30.2861 % austin 78712
% -95.6668	29.8297 % houston zip code 77084 
% -96.6745	32.7100 % dallas 75217 
% -98.525262	29.59002 % San antonip   78248

gzy=city(cityvalue,2); 
gzx=city(cityvalue,1);
boxsize=0.15;
boxboundary=[gzx+boxsize/2,gzx-boxsize/2,gzy+boxsize/2,gzy-boxsize/2]
maxlon=boxboundary(1);
minlon=boxboundary(2);	
maxlat=boxboundary(3);	
minlat=boxboundary(4);
%N=4000000;
count=0;
b=find((xbig(:,1+4)<maxlon)&(xbig(:,1+4)>minlon)&(xbig(:,2+4)<maxlat)&(xbig(:,2+4)>minlat));
xallinbox=xbig(b,:);
bbridges=find(xallinbox(:,3)==1);
xbridgesinbox=xallinbox(bbridges,:);
clearvars  xbig xdbig xd0 xd1 xd2 xd3 xd4 xd5 xd6 xd7 xd8 xd9 xd10 xd11 xd12 xd13 xd14 xd15 isbridge istunnel
%xallinbox columes are 1: maxspeed(usually rubbish) 2: layer ( whatever
%this means, a lot of rubbish here) 3: isbridge 4: istunnel 5: start lon 6
%start lat 7: end lon 8: endlat
disp('Variables cleared');toc;
    
x=xallinbox(:,5:8);
xbridge=xbridgesinbox(:,5:8);
disp('road box done');toc;
disp(['Number of bridges=' num2str(length(xbridge(:,1))) ]);
%csvwrite('C:\Users\vaidy\Documents\Final_Road_data\Final_Road_data\pennsylvania\pitt_1.csv',x);
%extract the node from the edgelist data
xnodes=unique([x(:,1) x(:,2); x(:,3) x(:,4)],'rows');
xnodesbridges=unique([xbridge(:,1) xbridge(:,2); xbridge(:,3) xbridge(:,4)],'rows');
[xedges ix ixedges]=unique([x(:,1) x(:,2) x(:,3) x(:,4)],'rows');
xcounts=accumarray(ixedges,1);
edgelength=111*power(power((x(:,1)-x(:,3)),2)+power((x(:,2)-x(:,4)),2),0.5);
edgelength_bridges=111*power(power((xbridge(:,1)-xbridge(:,3)),2)+power((xbridge(:,2)-xbridge(:,4)),2),0.5);
%x=xbig(1:N,:);
%return
edgelabels=generatenodeid(x,xnodes);
%bridgelabels=generatenodeid(xbridge,xnodesbridges);
bridgelabels=edgelabels(bbridges);

xnodelabels=1:1:length(xnodes);
xnodelabels=xnodelabels(:);



%txcs= csvread('C:\Users\vaidy\Google Drive\TINA\us310mcc.csv'); 
%txptx=substationloc();
disp('dataread done');toc;

xwithlengths=[x edgelength];
xboundaryindex=boundary(xnodes(:,1),xnodes(:,2));
boundaryperimeter=111*perimeter(alphaShape(xnodes(:,1),xnodes(:,2)));
xboundary=[xnodes(xboundaryindex,1) xnodes(xboundaryindex,2)];

bitarrayforlocationsofbridges=ismember(xnodes,xnodesbridges,'rows');
bridgenodelocs=find(bitarrayforlocationsofbridges==1);


%return
%choose S nodes wisely
%snodes=xnodelabels(I2);
snodes=xnodelabels(xboundaryindex);
edgelength_base=edgelength/.0001;%normalized edgelength becasue in matal zero means disconnected
G_base=graph(edgelabels(:,1),edgelabels(:,2),edgelength_base,length(xnodelabels));
edgelength_nobridges=edgelength.*(1-xallinbox(:,3)); % zero means disconneted in MATLAB graph model
G=graph(edgelabels(:,1),edgelabels(:,2),edgelength_nobridges,length(xnodelabels));
%A_G=adjacency(G);
%P_G=A_G./(sum(A_G,2));

sps1_base=distances(G_base,xboundaryindex(1:end-1,:),'Method','positive');
sps2_base=distances(G_base,xboundaryindex(1:end-1,:),unique(snodes),'Method','positive');
%sps2_base(sps2_base==inf)  = max(sps2_base(isfinite(sps2_base)));
sps2rawmax_base=max(max(sps2_base(isfinite(sps2_base))));
sps2_base(sps2_base==inf)  = boundaryperimeter/2;



% 
 sps1=distances(G,xboundaryindex(1:end-1,:),'Method','positive');
 sps2=distances(G,xboundaryindex(1:end-1,:),unique(snodes),'Method','positive');
 %sps2(sps2==inf)  = max(sps2(isfinite(sps2)));
 sps2rawmax=max(max(sps2(isfinite(sps2))));
 sps2(sps2==inf)  = boundaryperimeter/2;

delta_edgelength=edgelength_base-edgelength_nobridges;
TI_base=min(sps2_base);
TA_base=mean(sps2_base);
TM_base=max(sps2_base);
TD_base=3*TA_base-(TI_base+TM_base);% likelihood from triangle model


TI=min(sps2);
TA=mean(sps2);
TM=max(sps2);
TD=3*TA-(TI+TM);% likelihood from triangle model
toc
delta_TI=TI-TI_base;
delta_TA=TA-TI_base;
delta_TM=TM-TI_base;
delta_TD=TD-TI_base;
toc
disp('delta done');
plotthedeltas(delta_TI,delta_TA,delta_TM, outfolder);
%return
betweennessofnodes_base=centrality(G_base,'betweenness'); toc
betweennessofnodes=centrality(G,'betweenness');toc



xnodebridgeslabels=xnodelabels(bridgenodelocs);

betweennessofbridgenodes=betweennessofnodes_base(bridgenodelocs);



delta_betweennessofnodes=betweennessofnodes-betweennessofnodes_base;
toc

%bins = conncomp(G);

xnodeswitbridgeflag=[xnodes bitarrayforlocationsofbridges betweennessofnodes_base];
xnodeswitbridgeflagsorted=sortrows(xnodeswitbridgeflag,4,'descend');%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% modified the graph in sequential order of bridges and recalculate the
%%%% betweennesscentrality : THIS code is is complicated
% betweeness is for node. so how to remove edges???lol
% create h edge table with node betweness for start and end nodes and then
% the edge betweenes as a fucntion of the node betweenenss
[lia_S, locb_S]=ismember(xallinbox(:,5:6), xnodeswitbridgeflag(:,1:2),'rows');
[lia_E, locb_E]=ismember(xallinbox(:,7:8), xnodeswitbridgeflag(:,1:2),'rows');
b_S=xnodeswitbridgeflag(locb_S,end);
b_E=xnodeswitbridgeflag(locb_E,end);
b_SE=(b_S+b_E)/2;
b_deltaSE=(b_S-b_E);

xallinbox=[xallinbox b_S b_E b_SE b_deltaSE];

xallinboxwithbridgesorted=sortrows(xallinbox,[3 11 12],'descend');
N_bridge=length(xallinboxwithbridgesorted(find(xallinboxwithbridgesorted(:,3)==1),12));
xnodeswitbridgesorted=unique(xnodeswitbridgeflagsorted,'rows');
xboundaryindex=boundary(xnodeswitbridgesorted(:,1),xnodeswitbridgesorted(:,2));
snodes=xnodelabels(xboundaryindex);
%Tauparam=zeros(N_bridge,3);
numberofedges=length(xallinboxwithbridgesorted(:,1));
TI_ithbridgeremoved=zeros(N_bridge,length(unique(snodes)));
TA_ithbridgeremoved=zeros(N_bridge,length(unique(snodes)));
TM_ithbridgeremoved=zeros(N_bridge,length(unique(snodes)));
tic
%% Delta with all but k bridges removed
bridgekickfraction=0;% 0 means 1 bridge will be kicvked. becaseu no brdige being kicked is irrelvant in the plot 1 mean all bridges will be kicked
k=bridgekickfraction*N_bridge;
v_allbutone=zeros(numberofedges,1);
v_allbutone(1:N_bridge-k,:)=1;


G_allbutone=graph(edgelabels(:,1),edgelabels(:,2),edgelength,length(xnodelabels));
%sps4=distances(G_allbutone,xboundaryindex(1:end-1,:),unique(snodes),'Method','positive');
sps4base=distances(G_allbutone,unique(snodes),unique(snodes),'Method','positive');
sps4baserawmax=max(max(sps4base(isfinite(sps4base))));
sps4base(sps4base==inf)=boundaryperimeter/2;
TI_sps4base=min(sps4base);
TA_sps4base=mean(sps4base);
TM_sps4base=max(sps4base);

edgelength_withouthithbridge=edgelength./(1-v_allbutone);



G_allbutone=graph(edgelabels(:,1),edgelabels(:,2),edgelength_withouthithbridge,length(xnodelabels));
%sps4=distances(G_allbutone,xboundaryindex(1:end-1,:),unique(snodes),'Method','positive');
sps4=distances(G_allbutone,unique(snodes),unique(snodes),'Method','positive');
%sps4(sps4==inf)  = max(sps4(isfinite(sps4)));
sps4rawmax=max(max(sps4(isfinite(sps4))));
sps4(sps4==inf)=boundaryperimeter/2;
%T_i=minspantree(G_i);
TI_AllButOnebridgeremoved=min(sps4);
TA_AllButOnebridgeremoved=mean(sps4);
TM_AllButOnebridgeremoved=max(sps4);
delta_TI_AllButOnebridgeremoved=TI_AllButOnebridgeremoved-TI_sps4base;
delta_TA_AllButOnebridgeremoved=TA_AllButOnebridgeremoved-TA_sps4base;
delta_TM_AllButOnebridgeremoved=TM_AllButOnebridgeremoved-TM_sps4base;


figure(100000);
hold on;
scatter(1:1:length(TI_AllButOnebridgeremoved),TI_AllButOnebridgeremoved,'g','s');
scatter(1:1:length(TA_AllButOnebridgeremoved),TA_AllButOnebridgeremoved,'b','o');
scatter(1:1:length(TM_AllButOnebridgeremoved),TM_AllButOnebridgeremoved,'r','p');

scatter(1:1:length(TM_sps4base),TI_sps4base,'g','s','filled');
scatter(1:1:length(TM_sps4base),TA_sps4base,'b','o','filled');
scatter(1:1:length(TM_sps4base),TM_sps4base,'r','p','filled');

xlabel('boundary node id');
ylabel('change in delay parameters');
legend('Min  bridge removed','Mean  bridge removed','Max bridge removed','Min','Mean','Max','Location', 'southeast');
title(['BaseandKRemovedDelays k=' num2str(k+1) 'cityval= ' num2str(cityvalue)  ]);
print(strcat(outfolder, 'BaseandKRemovedDelays'), "-dpng");
savefig(strcat(outfolder, 'BaseandKRemovedDelays'));
hold off
toc
%return

figure(100001);
hold on;
scatter(1:1:length(delta_TI_AllButOnebridgeremoved),delta_TI_AllButOnebridgeremoved,'g','s');
scatter(1:1:length(delta_TA_AllButOnebridgeremoved),delta_TA_AllButOnebridgeremoved,'b','o');
scatter(1:1:length(delta_TM_AllButOnebridgeremoved),delta_TM_AllButOnebridgeremoved,'r','p');
xlabel('boundary node id');
ylabel('change in delay parameters');
legend('Min','Mean','Max','Location', 'southeast');
title([ 'DelayChangeKBridges k=' num2str(k+1) 'cityval= ' num2str(cityvalue)  ]);
print(strcat(outfolder, 'DelayChangeKBridges'), "-dpng");
savefig(strcat(outfolder, 'DelayChangeKBridges'));
hold off
toc
return



%% ith bridge removed

v=zeros(numberofedges,1);
numberofbridgestoremove=N_bridge/100;
for i=1:1:numberofbridgestoremove
    %v=zeros(numberofedges,1);
    v(i)=1;
    edgelength_withouthithbridge=edgelength./(1-v);
    G_i=graph(edgelabels(:,1),edgelabels(:,2),edgelength_withouthithbridge,length(xnodelabels));
    sps3=distances(G_i,xboundaryindex(1:end-1,:),unique(snodes),'Method','positive');
    sps3(sps3==inf)  = max(sps3(isfinite(sps3)));
    %T_i=minspantree(G_i);
    TI_ithbridgeremoved(i,:)=min(sps3);
    TA_ithbridgeremoved(i,:)=mean(sps3);
    TM_ithbridgeremoved(i,:)=max(sps3);
    %Tauparam(i, :)=[TI_ithbirdgeremovesd TA_ithbirdgeremovesd TM_ithbirdgeremovesd];
    
    clearvars G_i 
    toc
end

%%
disp(['Number of nodes = ' num2str(length(xnodes(:,1)))]);
disp(['Number of edges = ' num2str(length(x(:,1)))]);
disp(['Number of bridges = ' num2str(length(xbridge(:,1)))]);
disp(['perimeter=' num2str(boundaryperimeter)]);
disp(['maximum of finite paths with '   num2str(k+1)  ' bridges removed= ' num2str(sps4rawmax)]);
disp(['maximum of finite paths with base road network  =' num2str(sps2rawmax)]);
%plotthedeltas(delta_TI,delta_TA,delta_TM, outfolder);
toc
return
Delta_TI_ithbridgeremoved=TI_ithbridgeremoved-TI_base;
Delta_TA_ithbridgeremoved=TA_ithbridgeremoved-TA_base;
Delta_TM_ithbridgeremoved=TM_ithbridgeremoved-TM_base;


csvwrite(strcat(outfolder, 'TI_ithbridgeremoved.csv'),TI_ithbridgeremoved);
csvwrite(strcat(outfolder, 'TA_ithbridgeremoved.csv'),TA_ithbridgeremoved);
csvwrite(strcat(outfolder, 'TM_ithbridgeremoved.csv'),TM_ithbridgeremoved);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;
%return
makevideo=0;

minxlim=min(xnodeswitbridgeflagsorted(:,1));
maxxlim=max(xnodeswitbridgeflagsorted(:,1));
minylim=min(xnodeswitbridgeflagsorted(:,2));
maxylim=max(xnodeswitbridgeflagsorted(:,2));
minclim=min(xnodeswitbridgeflagsorted(:,4));
maxclim=max(xnodeswitbridgeflagsorted(:,4));
switch makevideo 
   case 1
    topN=5000;
    videodir=strcat("C:\Users\vaidy\Documents\bridgevideos\", thisrun,"\");
    mkdir(videodir)
    for i=1:length(xnodeswitbridgeflagsorted(1:topN,1))
        bridgesbuilt=sum(xnodeswitbridgeflagsorted(1:i,3));
        scatter(xnodeswitbridgeflagsorted(1:i,1),xnodeswitbridgeflagsorted(1:i,2),10,xnodeswitbridgeflagsorted(1:i,4));colormap(jet);colorbar;
        title(['node= ' num2str(i) ' bridge node= ' num2str(bridgesbuilt)]);
        set(gca,'XLim',[minxlim, maxxlim],'YLim',[minylim, maxylim],'CLim',[minclim, maxclim]);
        print(strcat(videodir, num2str(i)), "-dpng")
    end
    otherwise
end
%return
%% Generate all the plots
%scatter(xnodeswitbridgeflagsorted(1:1400,1),xnodeswitbridgeflagsorted(1:1400,2));
%plotthedeltas(delta_TI,delta_TA,delta_TM, thisrun);
plotthedeltas(delta_TI,delta_TA,delta_TM, outfolder);
return;
%betweennessplots(betweennessofnodes_base,xnodes,strcat("base" , thisrun));
betweennessplots(betweennessofnodes_base,xnodes,strcat(outfolder,"\","base" ));
%betweennessplots(betweennessofnodes,xnodes, thisrun);
betweennessplots(betweennessofnodes,xnodes, outfolder);
%betweennessplots(betweennessofbridgenodes,xnodesbridges, strcat("only bridges" , thisrun))
betweennessplots(betweennessofbridgenodes,xnodesbridges, strcat(outfolder,"\","only bridges"))

%betweennessplots(delta_betweennessofnodes,xnodes, thisrun);toc
betweennessplots(delta_betweennessofnodes,xnodes, outfolder);toc

plotbridgereusltsall(strcat(outfolder,"\","base"),edgelength_base,edgelength_bridges,snodes,xnodes,TI_base,TD_base,TM_base,TA_base,sps2_base,xboundary,xbridge,x);
plotbridgereusltsall(outfolder,edgelength,edgelength_bridges,snodes,xnodes,TI,TD,TM,TA,sps2,xboundary,xbridge,x);

%plotbridgereusltsall(strcat("base" , thisrun),edgelength_base,edgelength_bridges,snodes,xnodes,TI_base,TD_base,TM_base,TA_base,sps2_base,xboundary,xbridge,x);
%plotbridgereusltsall(thisrun,edgelength,edgelength_bridges,snodes,xnodes,TI,TD,TM,TA,sps2,xboundary,xbridge,x);


toc

%%

function betweennessplots(betweennessofbridgenodes,xnodes, outfolder)
    figure();
    histogram(log10(betweennessofbridgenodes));
    xlabel('log10(nodebetweenness)');
    ylabel('counts');    
    print(strcat(outfolder, 'histofbetweenness'), "-dpng");
    figure();
    scatter(xnodes(:,1),xnodes(:,2),1,log10(betweennessofbridgenodes),'*');
    colormap(jet)
    colorbar
    xlabel('longitude');
    ylabel('latitude');
    title(num2str(toc));
    print(strcat(outfolder, 'betweennesssmap'), "-dpng");
    savefig(strcat(outfolder, 'betweennesssmap'));
end
function plotthedeltas(delta_TI,delta_TA,delta_TM, outfolder)
    hold off;
figure();
hold on;
% stem(delta_TI);
% stem(delta_TA);
% stem(delta_TM);
scatter(1:1:length(delta_TI),delta_TI,'g','s');
scatter(1:1:length(delta_TA),delta_TA,'b','o');
scatter(1:1:length(delta_TM),delta_TM,'r','*');

xlabel('boundary node id');
ylabel('increase in delay');
legend('Min','Mean','Max','Location', 'southeast');
title('DelayChangeduetoNoBridges')
print(strcat(outfolder, 'DelayChangeduetoNoBridges'), "-dpng");
savefig(strcat(outfolder, 'DelayChangeduetoNoBridges'));
hold off
end

function plotthedeltas_somebridges(delta_TI,delta_TA,delta_TM, outfolder)
    hold off;
figure();
hold on;
% stem(delta_TI);
% stem(delta_TA);
% stem(delta_TM);
scatter(1:1:length(delta_TI),delta_TI,'g','s');
scatter(1:1:length(delta_TA),delta_TA,'b','o');
scatter(1:1:length(delta_TM),delta_TM,'r','*');

xlabel('boundary node id');
ylabel('increase in delay');
legend('Min','Mean','Max','Location', 'southeast');
title(num2str(toc))
print(strcat(outfolder, 'DelayChangeduetoIBridgesremoved'), "-dpng");
savefig(strcat(outfolder, 'DelayChangeduetoIBridgesremoved'));
hold off
end

function plotbridgereusltsall(thisrun,edgelength,edgelength_bridges,snodes,xnodes,TI,TD,TM,TA,sps2,xboundary,xbridge,x)
    hold on;
    figure();
    fignum=0;
    xnodesize=1;

    cnodesize=50;
    pnodesize=150;
    scatter(xnodes(:,1), xnodes(:,2),'.');
    xlabel('longitude');
    ylabel('latitude');
    print(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetnodes", thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetnodes", thisrun));
    
    print(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetnodes", thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetnodes", thisrun));

    fignum=fignum+1;figure(fignum)
    plot(xboundary(:,1),xboundary(:,2),'red','LineWidth',2);
    scatter(xboundary(:,1),xboundary(:,2),'black','*','LineWidth',2)
    Xs=x(:,1); Xe=x(:,3); Ys=x(:,2); Ye=x(:,4); X= [Xs Xe]; Y=[Ys Ye];
    Xsbridge=xbridge(:,1); Xebridge=xbridge(:,3); Ysbridge=xbridge(:,2); Yebridge=xbridge(:,4); Xbridge= [Xsbridge Xebridge]; Ybridge=[Ysbridge Yebridge];
    line(X',Y','color','black','LineWidth',.1);
    %legend('road');
    line(Xbridge',Ybridge','color','green','LineWidth',.1);
    %legend('bridge');
    xlabel('longitude');
    ylabel('latitude');
    title(num2str(toc))
    print(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednet", thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednet", thisrun));

    fignum=fignum+1;figure(fignum)
    plot(xboundary(:,1),xboundary(:,2),'red','LineWidth',2);
    scatter(xboundary(:,1),xboundary(:,2),'black','*','LineWidth',2)
    %Xs=x(:,1); Xe=x(:,3); Ys=x(:,2); Ye=x(:,4); X= [Xs Xe]; Y=[Ys Ye];
    Xsbridge=xbridge(:,1); Xebridge=xbridge(:,3); Ysbridge=xbridge(:,2); Yebridge=xbridge(:,4); Xbridge= [Xsbridge Xebridge]; Ybridge=[Ysbridge Yebridge];
    %line(X',Y','color','black','LineWidth',.1);
    %legend('road');
    line(Xbridge',Ybridge','color','green','LineWidth',.1);
    %legend('bridge');
    xlabel('longitude');
    ylabel('latitude');
    title(num2str(toc))
    print(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetbridgesonly", thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetbridgesonly", thisrun));
    
    fignum=fignum+1;figure(fignum)
    hold on
    plot(xboundary(:,1),xboundary(:,2),'red','LineWidth',2)
    scatter(xnodes(:,1),xnodes(:,2),xnodesize,'black','.');
    %scatter(ccs(:,1),ccs(:,2), cnodesize,'red','s');
    %scatter(ptx(:,1),ptx(:,2), cnodesize,'blue','*'); 
    xlabel('longitude');
    ylabel('latitude');
    title(num2str(toc))
    print(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetnodes", thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\undamagednetnodes", thisrun));
    
    fignum=fignum+1;figure(fignum)
    hold on
    histogram(1000*edgelength);
    histogram(1000*edgelength_bridges);
    xlabel('edgelength(meters)');
    ylabel('frequency');
    title('edgelegth distribution')
    title(num2str(toc))
    legend('all', 'bridges');
    %savefig(['C:\Users\vaidy\bridges\' 'yieldMT' num2str(Y)  'fullnetworkdamaged' '.fig'],'compact');
    print(strcat("C:\Users\vaidy\Google Drive\bridges\edgelength" , thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\edgelength" , thisrun));

    fignum=fignum+1;figure(fignum)
    hold on
    histogram(log10(1000*edgelength),'facealpha',.75); hold on;
    histogram(log10(1000*edgelength_bridges),'facealpha',.75); hold on;
    xlabel('log_{10}edgelength(meters)');
    ylabel('frequency');
    title(['edgelength distribution', num2str(toc)])
    legend('all', 'bridges');
    %savefig(['C:\Users\vaidy\bridges\' 'yieldMT' num2str(Y)  'fullnetworkdamaged' '.fig'],'compact');
    print(strcat("C:\Users\vaidy\Google Drive\bridges\edgelengthlog10" , thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\edgelengthlog10" , thisrun));


    forplotxsnodes=unique(snodes);
    fignum=fignum+1;figure(fignum)
    hold on
    %line(X',Y','color','black','LineWidth',0.1);     
    scatter(xnodes(forplotxsnodes,1),xnodes(forplotxsnodes,2),cnodesize,TI,'*');
    xlabel('longitude');
    ylabel('latitude');
    title('Min time T_i');
    colormap(jet)
    colorbar
    set(gca,'XLim',[min(xnodes(:,1)), max(xnodes(:,1))],'YLim',[min(xnodes(:,2)), max(xnodes(:,2))],'CLim',[min(min(sps2)), max(max(sps2))]);
    %savefig(['C:\Users\vaidy\bridges\' 'yieldMT' num2str(Y)  'fullnetworkdamaged' '.fig'],'compact');
    print(strcat("C:\Users\vaidy\Google Drive\bridges\Ti" , thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\Ti" , thisrun));
    
    fignum=fignum+1;figure(fignum)
    hold on
    %line(X',Y','color','black','LineWidth',0.05);     

    scatter(xnodes(forplotxsnodes,1),xnodes(forplotxsnodes,2),cnodesize,TD,'*');

    xlabel('longitude');
    ylabel('latitude');
    title('Most likely time T_d');
    colormap(jet)
    colorbar
    set(gca,'XLim',[min(xnodes(:,1)), max(xnodes(:,1))],'YLim',[min(xnodes(:,2)), max(xnodes(:,2))],'CLim',[min(min(sps2)), max(max(sps2))]);
    %savefig(['C:\Users\vaidy\bridges\' 'yieldMT' num2str(Y)  'fullnetworkdamaged' '.fig'],'compact');
    print(strcat("C:\Users\vaidy\Google Drive\bridges\Td" , thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\Td" , thisrun));

    fignum=fignum+1;figure(fignum)
    hold on
    %line(X',Y','color','black','LineWidth',0.05);     
    scatter(xnodes(forplotxsnodes,1),xnodes(forplotxsnodes,2),cnodesize,TM,'*');
    xlabel('longitude');
    ylabel('latitude');
    title('Max time T_m');
    colormap(jet)
    colorbar
    set(gca,'XLim',[min(xnodes(:,1)), max(xnodes(:,1))],'YLim',[min(xnodes(:,2)), max(xnodes(:,2))],'CLim',[min(min(sps2)), max(max(sps2))]);
    %savefig(['C:\Users\vaidy\bridges\' 'yieldMT' num2str(Y)  'fullnetworkdamaged' '.fig'],'compact');
    print(strcat("C:\Users\vaidy\Google Drive\bridges\Tm" , thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\Tm" , thisrun));

    fignum=fignum+1;figure(fignum)
    hold on
    scatter(xnodes(forplotxsnodes,1),xnodes(forplotxsnodes,2),cnodesize,TA,'*');
    %line(X',Y','color','black','LineWidth',0.05);     
    xlabel('longitude');
    ylabel('latitude');
    title('Mean time T_\mu');
    colormap(jet)
    colorbar
    set(gca,'XLim',[min(xnodes(:,1)), max(xnodes(:,1))],'YLim',[min(xnodes(:,2)), max(xnodes(:,2))],'CLim',[min(min(sps2)), max(max(sps2))]);
    %savefig(['C:\Users\vaidy\bridges\' 'yieldMT' num2str(Y)  'fullnetworkdamaged' '.fig'],'compact');
    print(strcat("C:\Users\vaidy\Google Drive\bridges\Ta" , thisrun), "-dpng");
    savefig(strcat("C:\Users\vaidy\Google Drive\bridges\Ta" , thisrun));
end







