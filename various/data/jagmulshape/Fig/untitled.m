clear; clc;close all
clear;clc;close all;
%% Plot the case of wider shot coverage       
zp=[0 10] ;
vp=[4100   4100];%velocity model
zp=0:10:2000;
vp=1800+3.1*zp;
zsrc=50;zrec=250;zd=2500;%source receiver and reflector depths;%source receiver and reflector depths
caprad=30;itermax=4;%offsets, cap radius, and max iter
pfan=-1;optflag=1;pflag=1;dflag=2;%default ray fan, and various flags

xstart1=-300;
xstart2=3200;

xrec=700:100:5200; 
color1='k';color2='k';color3='k';color4='k';

% plot receivers
% line(xrec,zrec*ones(size(xrec)),'color','k','linewidth',2,'linestyle','none','marker','v')
% 
% plot sources

for k =1:12
    for i=k+5:-5:k
        line(-i,k,'color','k','linewidth',5,'linestyle','none','marker','.','markersize',35)
    end
end
annotation('line',[0.07,0.1],[0.98 0.98],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
annotation('line',[0.07,0.07],[0.08 0.981],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
annotation('line',[0.07,0.1],[0.08 0.08],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
annotation('line',[0.92,0.95],[0.98 0.98],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
annotation('line',[0.95,0.95],[0.08 0.981],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
annotation('line',[0.92,0.95],[0.08 0.08],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
% annotation('line',[0.02,0.05],[0.98 0.98],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
% annotation('line',[0.022,0.022],[0.03 0.981],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
% annotation('line',[0.02,0.05],[0.03 0.03],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
% annotation('line',[0.96,0.99],[0.98 0.98],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
% annotation('line',[0.99,0.99],[0.03 0.981],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]
% annotation('line',[0.96,0.99],[0.03 0.03],'Linewidth',2,'color','k') %[begx,endx],[begy,endy]

set(gca,'fontname','times new roman');

set(gcf,'Units','centimeters','Position',[11 5 15 15]);
set(gca,'fontsize',20);
set(get(gca,'xlabel'),'fontsize',20)
set(get(gca,'ylabel'),'fontsize',20)
axis([-11.5 -6.5  0.5 11.5])
axis off
set(gcf,'Units','centimeters','Position',[11 5 15 15]);
set(gca,'Position',[0.01 0.01 0.99 0.99]);
% text(4000,2000,'Denser shot coverage system','fontsize',20)