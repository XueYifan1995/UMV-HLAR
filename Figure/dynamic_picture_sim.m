clc,clear,close all

%% Importing data
load npsAUV_zigzag_1010_005 npsAUV_zigzag_1010_005
load npsAUV_zigzag_2020_005 npsAUV_zigzag_2020_005
load npsAUV_zigzag_2505_005 npsAUV_zigzag_2505_005
load npsAUV_zigzag_2510_005 npsAUV_zigzag_2510_005

xinput = [npsAUV_zigzag_1010_005,npsAUV_zigzag_2020_005];
x_1_self = xinput(1:3,:);
u = xinput(4,:)*pi/180;


%% Plotting raw data
figure
plot(xinput(1,:),'linewidth',6,'color',[0,0.45,0.74])  %u
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim1u','-r600','-dpng');  

figure
plot(xinput(2,:),'linewidth',6,'color',[0.93,0.69,0.13])  %v
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim2v','-r600','-dpng');  %print figures

figure
plot(xinput(3,:),'linewidth',6,'color',[0.47,0.67,0.19])  %r
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim3r','-r600','-dpng');  

figure
plot(xinput(4,:),'linewidth',6,'color',[0,0,0])  %d
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim4d','-r600','-dpng');  

figure
plot(xinput(1,:).*abs(xinput(1,:)),'linewidth',6,'color',[0,0.45,0.74])  %u|u|
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim5uabsu','-r600','-dpng');  

figure
plot(xinput(1,:).*xinput(2,:),'linewidth',6,'color',[1.00,0.41,0.16])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim6uv','-r600','-dpng');  

figure
plot(xinput(1,:).*xinput(3,:),'linewidth',6,'color',[1.00,0.41,0.16])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim7ur','-r600','-dpng');  

figure
plot(xinput(2,:).*xinput(2,:),'linewidth',6,'color',[84/255,130/255,53/255])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim8vv','-r600','-dpng');  

figure
plot(xinput(2,:).*xinput(3,:),'linewidth',6,'color',[84/255,130/255,53/255])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim9vr','-r600','-dpng');  

figure
plot(xinput(2,:).*abs(xinput(2,:)),'linewidth',6,'color',[1.00,0.41,0.16])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim10vabsv','-r600','-dpng');  

figure
plot(xinput(3,:).*xinput(3,:),'linewidth',6,'color',[84/255,130/255,53/255])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim11rr','-r600','-dpng');  

figure
plot(xinput(3,:).*abs(xinput(3,:)),'linewidth',6,'color',[1.00,0.41,0.16])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim12rabsr','-r600','-dpng');  

figure
plot(xinput(1,:).*xinput(1,:).*xinput(4,:),'linewidth',6,'color',[1.00,0.41,0.16])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim13uud','-r600','-dpng'); 

figure
plot(xinput(3,:).*xinput(3,:).*xinput(3,:),'linewidth',6,'color',[1.00,0.41,0.16])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim14rrr','-r600','-dpng');  

figure
plot(xinput(1,:).*xinput(1,:).*xinput(4,:).*xinput(4,:),'linewidth',6,'color',[84/255,130/255,53/255])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim15uudd','-r600','-dpng');  

figure
plot(1./xinput(1,:),'linewidth',6,'color',[0,0.45,0.74])  %uv
box off
axis off
set(gcf,'unit','normalized','position',[0.2,0.2,1,0.12]);
print('Sim161u','-r600','-dpng');  