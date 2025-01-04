%% FYI coring, refined temperatures
close all; clear; clc; 
load('T66.mat',"t_T66","T_T66","lat_T66","long_T66"); % Thermistor 2019T66 data import (FYI)
cut = 1097; t_T66 = t_T66(1:cut); T_T66 = T_T66(1:cut,:); lat_T66 = lat_T66(1:cut); lon_T66 = long_T66(1:cut);
td_0_T66 = (datenum(t_T66)- datenum(t_T66(1))); td_T66 = datenum(t_T66); % time of experiment [d]
top_t_dir_T66 = [0 212 215 223 225 233 241 246 254 262 266 269 272]; % Direct time of measurements [d]
top_dir_T66 =  -[0   0   6  12   0   2   0  12  22  40  42  52  54]/100; % Thickness for interpolation [m]
bot_t_dir_T66 = [ 0  4  5 24 45 64  99 154 195 212 244 261 272]; % Direct time of measurements [d]
bot_dir_T66 =  -[33 40 42 54 66 82 104 140 154 154 146 138 126]/100; % Thickness for interpolation [m]
bot_int_T66 = interp1(bot_t_dir_T66,bot_dir_T66,td_0_T66,'pchip'); % Experimental bottom interpolation
top_int_T66 = interp1(top_t_dir_T66,top_dir_T66,td_0_T66,'pchip'); % Experimental bottom interpolation
z_T66 = 0.96:-0.02:-4.8+0.96; % IMB vertical coordinates, original
z_T66_ref = 0.0:-0.005:-1.8; % IMB vertical coordinates, refined, 5 mm spacing
T_T66_ref = zeros(size(T_T66,1),length(z_T66_ref));
for i = 1:length(T_T66)
    T_T66_ref(i,:) = interp1(z_T66,T_T66(i,:),z_T66_ref,'linear'); % refinement of IMB temperature
end
frz_T66_ref = T_T66_ref; % temperature of only frozen sesnors
for i = 1:size(T_T66_ref,1) % time
    for j = 1:size(T_T66_ref,2) % depth
        if  j > (z_T66_ref(1) - bot_int_T66(i))/-diff(z_T66_ref(1:2))
            frz_T66_ref(i,j) = NaN;
        elseif  j < (z_T66_ref(1) - top_int_T66(i))/-diff(z_T66_ref(1:2))
            frz_T66_ref(i,j) = NaN;
        end
    end
end
clearvars cut td_0_T66 top_t_dir_T66 top_dir_T66 bot_t_dir_T66 bot_dir_T66 T_T66 i j long_T66

% figure
% range = [-30 -25 -20 -15 -10 -8 -6 -4 -2 0 2 4]; 
% contourf(datenum(t_T66),z_T66_ref,frz_T66_ref',range,'-','ShowText','on','LabelSpacing',400,'LineColor','none','LineWidth',0.1); hold on
% plot(datenum(t_T66),top_int_T66,'k--'); plot(datenum(t_T66),bot_int_T66,'k--');
% hYLabel = ylabel('Depth (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
% clim([-30 5]);
% yticks(-1.8:0.2:0.2); ylim([-1.8 0.2]); % depth limits
% colorbar; brighten(.1); colormap(gca,bluewhitered);
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% hBar1 = colorbar; ylabel(hBar1,'Ice temperature (°C)','FontSize',8);
% clearvars range hYLabel ax hBar1

%% Stretch of T from IMB towards surface interface
T_fy_int = frz_T66_ref'; zT_fy_int = z_T66_ref;
for i = 1:length(T_fy_int)
     zT_fy_top(i) = zT_fy_int(find(~isnan(T_fy_int(:,i)),1));
     zT_fy_bot(i) = zT_fy_int(find(~isnan(T_fy_int(:,i)),1,'last'));
     zT_fy_top_sur(i) = zT_fy_top(i)-zT_fy_top(i);
     zT_fy_bot_sur(i) = zT_fy_bot(i)-zT_fy_top(i);
end
zT_fy_str = zT_fy_int(:)*ones(1,length(T_fy_int));
for i = 1:length(T_fy_int)
    zT_fy_str(:,i) = (zT_fy_str(:,i)-zT_fy_top(i))*(+zT_fy_bot_sur(i)-zT_fy_top_sur(i))/(zT_fy_bot(i)-zT_fy_top(i)) + zT_fy_top_sur(i);
end
T_fy_int_sur = T_fy_int; % preallocation
for i = 1:length(T_fy_int)
    T_fy_int_sur(:,i) = interp1(zT_fy_str(:,i),T_fy_int(:,i),zT_fy_int','linear'); % interpolation of IMB temp along salinity core depth
end
[X,Y] = meshgrid(datenum(t_T66),zT_fy_int); Z = T_fy_int_sur; % contour plot preparation
clearvars i

figure
range = [-30 -25 -20 -15 -10 -8 -6 -4 -2 0 2 4]; 
contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','none','LineWidth',0.1); hold on
plot(datenum(t_T66),zT_fy_top-zT_fy_top,'k--'); plot(datenum(t_T66),zT_fy_bot-zT_fy_top,'k--');
hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
clim([-30 5]);
yticks(-1.8:0.2:0.2); ylim([-1.8 0.2]); % depth limits
colorbar; brighten(.1); colormap(gca,bluewhitered);
ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
hBar1 = colorbar; ylabel(hBar1,'Ice temperature (°C)','FontSize',8);
clearvars X Y Z range hYLabel ax hBar1

%% Salinity base data
load('Coring_old_fb.mat',"S_fy","zS_fy","zzS_fy","t_fy","td_fy","td"); % Coring data import
for i = 1:length(t_fy)
    [~,idx(i)] = min(abs(datenum(t_fy(i))-datenum(t_T66)));
    l_imb(i) = -zT_fy_bot(idx(i))+zT_fy_top(idx(i)); % core length from IMB for coring dates
end
lvl = [1:3 5:12 15 17:23]; % removing rafted FYI
dS_fy = zS_fy;
for i = 1:length(S_fy)
    hS(i) = zzS_fy{i,1}(end,2); % salinity core length
    dS_fy{i,1} = -zS_fy{i,1}*l_imb(i)/hS(i); % length of salinity measurements for FYI
end
hS(22) = zzS_fy{22,1}(end-3,2); dS_fy{22,1} = -zS_fy{22,1}*l_imb(22)/hS(22); % remove fb
hS(21) = zzS_fy{21,1}(end-2,2); dS_fy{21,1} = -zS_fy{21,1}*l_imb(21)/hS(21); % remove fb
hS(20) = zzS_fy{20,1}(end-1,2); dS_fy{20,1} = -zS_fy{20,1}*l_imb(20)/hS(20); % remove fb

S_lvl = S_fy; S_lvl(16)=[]; S_lvl(14)=[]; S_lvl(13)=[]; S_lvl(4)=[];
dS_fy_lvl = dS_fy; dS_fy_lvl(16)=[]; dS_fy_lvl(14)=[]; dS_fy_lvl(13)=[]; dS_fy_lvl(4)=[]; 
td_lvl = td_fy; td_lvl(16)=[]; td_lvl(14)=[]; td_lvl(13)=[]; td_lvl(4)=[];

S_lvl{18}(28:30) = []; dS_fy_lvl{18}(28:30) = []; td_lvl{18}(28:30) = []; % remove fb
S_lvl{17}(29:30) = []; dS_fy_lvl{17}(29:30) = []; td_lvl{17}(29:30) = []; % remove fb
S_lvl{16}(32) = []; dS_fy_lvl{16}(32) = []; td_lvl{16}(32) = [];  % remove fb

S_vect = cell2mat(S_lvl); dS_fy_vect = cell2mat(dS_fy_lvl); tdS_fy_vect = cell2mat(td_lvl);
x = tdS_fy_vect; y = dS_fy_vect; z = S_vect;
xv = datenum(t_T66);
yv = z_T66_ref;
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);
S_fy_int = Z;
clearvars l_imb lvl tdS_fy_vect dS_fy_vect S_vect idx x y z xv yv i S_lvl td td_fy t_fy hS td_lvl dS_fy_lvl zzS_fy

% figure
% range = [0 1 2 3 4 6 8 10]; % Salinity graph accuracy;
% [C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
% plot(datenum(t_T66),zT_fy_top-zT_fy_top,'k--'); plot(datenum(t_T66),zT_fy_bot-zT_fy_top,'k--');
% clabel(C,h,'FontSize',6,'Color','k');
% hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); colorbar; colormap summer; brighten(.7);
% clim([0 8]);
% yticks(-1.8:0.2:0.2); ylim([-1.8 0.2]); % depth limits
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% hBar1 = colorbar; ylabel(hBar1,'Ice salinity','FontSize',8);
% clearvars X Y Z C h range hYLabel ax hBar1

%% Salinity stretch over IMB depth
for i = 1:length(S_fy_int)
     zS_fy_top(i) = z_T66_ref(find(~isnan(S_fy_int(:,i)),1));
     zS_fy_bot(i) = z_T66_ref(find(~isnan(S_fy_int(:,i)),1,'last'));
     zT_fy_top(i) = z_T66_ref(find(~isnan(T_fy_int_sur(:,i)),1));
     zT_fy_bot(i) = z_T66_ref(find(~isnan(T_fy_int_sur(:,i)),1,'last'));
end
zT_fy_str = z_T66_ref(:)*ones(1,length(T_fy_int));
for i = 1:length(T_fy_int)
    zT_fy_str(:,i) = (zT_fy_str(:,i)-zS_fy_top(i))*(+zT_fy_bot(i)-zT_fy_top(i))/(zS_fy_bot(i)-zS_fy_top(i)) + zT_fy_top(i);
end

S_fy_int_sur = T_fy_int; % preallocation
for i = 1:length(T_fy_int)
    S_fy_int_sur(:,i) = interp1(zT_fy_str(:,i),S_fy_int(:,i),z_T66_ref','linear'); % interpolation of salinity along IMB depth
end
[X,Y] = meshgrid(datenum(t_T66),z_T66_ref); Z = S_fy_int_sur; % contour plot preparation
clearvars zT_fy_str i

figure
range = [0 1 2 3 4 6 8 10]; % Salinity graph accuracy;
[C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
plot(datenum(t_T66),zT_fy_top-zT_fy_top,'k--'); plot(datenum(t_T66),zT_fy_bot-zT_fy_top,'k--');
clabel(C,h,'FontSize',6,'Color','k');
hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
clim([0 8]);
yticks(-1.8:0.2:0.2); ylim([-1.8 0.2]); % depth limits
colorbar; colormap summer; brighten(.7);
ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
hBar1 = colorbar; ylabel(hBar1,'FYI salinity','FontSize',8);
clearvars X Y Z C range h hYLabel ax hBar1

%% Brine volume
vb_fy_int = zeros(size(S_fy_int_sur)); F1_fy_int = vb_fy_int;  F2_fy_int = vb_fy_int; rhoi_fy_int = vb_fy_int;
for j = 1:length(t_T66)
    for i = 1:length(z_T66_ref)
        if T_fy_int_sur(i,j) <= -2
            F1_fy_int(i,j) = -4.732-22.45*T_fy_int_sur(i,j) - 0.6397*T_fy_int_sur(i,j).^2 - 0.01074*T_fy_int_sur(i,j).^3; % Cox and Weeks (1983)
            F2_fy_int(i,j) = 8.903*10^-2 - 1.763*10^-2*T_fy_int_sur(i,j) - 5.33*10^-4*T_fy_int_sur(i,j)^2 - 8.801*10^-6*T_fy_int_sur(i,j)^3; % Cox and Weeks (1983)
        else
            F1_fy_int(i,j) = -4.1221*10^-2 + -1.8407*10^1*T_fy_int_sur(i,j) + 5.8402*10^-1*T_fy_int_sur(i,j)^2 + 2.1454*10^-1*T_fy_int_sur(i,j)^3; % F1 from Lepparanta and Manninen (1988)
            F2_fy_int(i,j) = 9.0312*10^-2 + -1.6111*10^-2*T_fy_int_sur(i,j) + 1.2291*10^-4*T_fy_int_sur(i,j)^2 + 1.3603*10^-4*T_fy_int_sur(i,j)^3; % F2 from Lepparanta and Manninen (1988)
        end
        rhoi_fy_int(i,j) = (917-1.403*10^-1*T_fy_int_sur(i,j)); % pure ice density, Pounder (1965) for T_insitu
        vb_fy_int(i,j) = (1000).*(rhoi_fy_int(i,j).*S_fy_int_sur(i,j)/1000)./(F1_fy_int(i,j)-rhoi_fy_int(i,j).*S_fy_int_sur(i,j)/1000.*F2_fy_int(i,j)); % Brine volume for T_insitu, LM
        if vb_fy_int(i,j) < 0
            vb_fy_int(i,j) = NaN;
        elseif vb_fy_int(i,j) > 1000
            vb_fy_int(i,j) = NaN;
        end
    end
end
clearvars i j F1_fy_int F2_fy_int rhoi_fy_int

% figure
% [X,Y] = meshgrid(datenum(t_T66),z_T66_ref); Z = vb_fy_int/1000; % contour plot preparation
% range = [0 0.02 0.05 0.1 0.2 0.3];
% contourf(X,Y,Z,range,'-','ShowText','off','LabelSpacing',400,'LineColor','none','LineWidth',0.1); hold on
% plot(datenum(t_T66),zT_fy_top-zT_fy_top,'k--'); plot(datenum(t_T66),zT_fy_bot-zT_fy_top,'k--');
% hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
% colorbar; colormap summer; brighten(0.5); clim([0 0.3]);
% text([td_T66(490) td_T66(800) td_T66(600) td_T66(900) td_T66(980) td_T66(980)],[-0.45 -1.38 -1.20 -0.9 -0.55 -0.30],{'0.02','0.1','0.05','0.1','0.2','0.3'},'color','k','FontSize',8);
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% yticks(-2.2:0.2:0.4); ylim([-1.8 0.2]); % summertime depth limits
% hBar1 = colorbar; ylabel(hBar1,'Brine volume','FontSize',8);
% clearvars X Y Z range hYLabel ax hBar1

%% Density base data
load('Coring_old_fb.mat',"rho_fy","zrho_fy","td_rho_fy","t_fy"); % Coring data import
for i = [1:3 5:6 8:23]
    [~,idx(i)] = min(abs(datenum(t_fy(i))-datenum(t_T66)));
    l_imb(i) = -zT_fy_bot(idx(i))+zT_fy_top(idx(i)); % core length from IMB for coring dates
end
for i = [1:3 5:6 8:23]
    h_rho(i) = zrho_fy{i,1}(end); % density core length
    zrho_fy{i,1} = -zrho_fy{i,1}*l_imb(i)/h_rho(i); % length of density measurements for FYI
end
rho_lvl = rho_fy; rho_lvl(7)=[]; rho_lvl(4)=[];
drho_fy_lvl = zrho_fy; drho_fy_lvl(7)=[]; drho_fy_lvl(4)=[]; 
td_lvl = td_rho_fy; td_lvl(7)=[]; td_lvl(4)=[];
drho_fy_lvl{19}(31:32)=[]; drho_fy_lvl{20}(30:31)=[]; rho_lvl{19}(31:32)=[]; rho_lvl{20}(30:31)=[]; td_lvl{19}(31:32)=[]; td_lvl{20}(30:31)=[]; % remove false bottom
drho_fy_lvl{21}(18:20)=[]; rho_lvl{21}(18:20)=[]; td_lvl{21}(18:20)=[]; % remove false bottom
drho_fy_lvl{1}(7:8)=[]; rho_lvl{1}(7:8)=[]; td_lvl{1}(7:8)=[]; % remove NaN week 1
drho_fy_lvl{2}(8:9)=[]; rho_lvl{2}(8:9)=[]; td_lvl{2}(8:9)=[]; % remove NaN week 2
drho_fy_lvl{3}(4:5)=[]; rho_lvl{3}(4:5)=[]; td_lvl{3}(4:5)=[]; % remove NaN week 3
drho_fy_lvl{4}(2)=[]; rho_lvl{4}(2)=[]; td_lvl{4}(2)=[]; % remove NaN week 4
rho_vect = cell2mat(rho_lvl);
drho_fy_vect = cell2mat(drho_fy_lvl);
tdvg_fy_vect = cell2mat(td_lvl);
x = tdvg_fy_vect; y = drho_fy_vect; z = rho_vect;
xv = datenum(t_T66); yv = z_T66_ref;
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y); rho_fy_int = Z;
clearvars l_imb lvl tdS_fy_vect dS_fy_vect S_vect idx x y z xv yv i S_lvl td td_fy t_fy hS td_lvl dS_fy_lvl zzS_fy

% figure
% range = [700 750 800 840 860 900 910 940]; % density graph accuracy
% [C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
% plot(datenum(t_T66),zT_fy_top-zT_fy_top,'k--'); plot(datenum(t_T66),zT_fy_bot-zT_fy_top,'k--');
% clabel(C,h,'FontSize',6,'Color','k');
% hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
% colorbar; colormap summer; brighten(0.5);
% clim([800 950]);
% yticks(-1.8:0.2:0.2); ylim([-1.8 0.2]); % depth limits
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% hBar1 = colorbar; ylabel(hBar1,'Ice density (kg/m^3)','FontSize',8);

%% Density stretch over IMB depth
for i = 1:length(rho_fy_int)
     zS_fy_top(i) = z_T66_ref(find(~isnan(rho_fy_int(:,i)),1));
     zS_fy_bot(i) = z_T66_ref(find(~isnan(rho_fy_int(:,i)),1,'last'));
     zT_fy_top(i) = z_T66_ref(find(~isnan(T_fy_int_sur(:,i)),1));
     zT_fy_bot(i) = z_T66_ref(find(~isnan(T_fy_int_sur(:,i)),1,'last'));
end
zT_fy_str = z_T66_ref(:)*ones(1,length(T_fy_int));
for i = 1:length(T_fy_int)
    zT_fy_str(:,i) = (zT_fy_str(:,i)-zS_fy_top(i))*(+zT_fy_bot(i)-zT_fy_top(i))/(zS_fy_bot(i)-zS_fy_top(i)) + zT_fy_top(i);
end

rho_fy_int_sur = T_fy_int; % preallocation
for i = 1:length(T_fy_int)
    rho_fy_int_sur(:,i) = interp1(zT_fy_str(:,i),rho_fy_int(:,i),z_T66_ref','linear'); % interpolation of salinity along IMB depth
end
[X,Y] = meshgrid(datenum(t_T66),z_T66_ref); Z = rho_fy_int_sur; % contour plot preparation
clearvars zT_fy_str i

figure
range = [700 750 800 840 860 900 910 940]; % density graph accuracy
[C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
plot(datenum(t_T66),zT_fy_top-zT_fy_top,'k--'); plot(datenum(t_T66),zT_fy_bot-zT_fy_top,'k--');
clabel(C,h,'FontSize',6,'Color','k');
hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
clim([800 950]);
yticks(-1.8:0.2:0.2); ylim([-1.8 0.2]); % depth limits
colorbar; colormap summer; brighten(0.5);
ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
hBar1 = colorbar; ylabel(hBar1,'Ice density (kg/m^3)','FontSize',8);
clearvars X Y Z C range h hYLabel ax hBar1

%% netCDF export
filename = 'MOSAiC_FYI_coring.nc';
time = datenum(t_T66-datetime('1979-01-01 00:00:00')); % time to nd datenum
delete(filename)

nccreate(filename,'time','Dimensions',{'time' length(time)},'format','netcdf4'); % time
ncwriteatt(filename,'time','standard_name','time');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','units','days since 1979-01-01 00:00:00');

nccreate(filename,'lat','Dimensions',{'time' length(time)});
ncwriteatt(filename,'lat','standard_name','lat');
ncwriteatt(filename,'lat','long_name','latitude');
ncwriteatt(filename,'lat','units','°');

nccreate(filename,'lon','Dimensions',{'time' length(time)});
ncwriteatt(filename,'lon','standard_name','lon');
ncwriteatt(filename,'lon','long_name','longitude');
ncwriteatt(filename,'lon','units','°');

nccreate(filename,'zi','Dimensions',{'zi' length(z_T66_ref)}); % Ice core length
ncwriteatt(filename,'zi','standard_name','zi');
ncwriteatt(filename,'zi','long_name','Ice core length');
ncwriteatt(filename,'zi','units','m');

nccreate(filename,'T_ice','Dimensions',{"zi",length(z_T66_ref),"time",length(time)}); % Temperature
ncwriteatt(filename,'T_ice','standard_name','Temperature');
ncwriteatt(filename,'T_ice','long_name','Ice temperature from 2019T66 buoy');
ncwriteatt(filename,'T_ice','units','°C');

nccreate(filename,'S_ice','Dimensions',{"zi",length(z_T66_ref),"time",length(time)}); % Salinity
ncwriteatt(filename,'S_ice','standard_name','Salinity');
ncwriteatt(filename,'S_ice','long_name','Ice practical salinity');
ncwriteatt(filename,'S_ice','units','PSU');

nccreate(filename,'vb','Dimensions',{"zi",length(z_T66_ref),"time",length(time)}); % Brine volume
ncwriteatt(filename,'vb','standard_name','Brine volume');
ncwriteatt(filename,'vb','long_name','Brine volume fraction');
ncwriteatt(filename,'vb','units','%');

nccreate(filename,'rho','Dimensions',{"zi",length(z_T66_ref),"time",length(time)}); % Ice density
ncwriteatt(filename,'rho','standard_name','Ice density');
ncwriteatt(filename,'rho','long_name','Ice density at -15°C laboratory temperature');
ncwriteatt(filename,'rho','units','kg/m3');

nccreate(filename,'hi','Dimensions',{'time' length(time)}); % Ice thickness
ncwriteatt(filename,'hi','standard_name','hi');
ncwriteatt(filename,'hi','long_name','Ice thickness from 2019T66 buoy');
ncwriteatt(filename,'hi','units','m');

ncwrite(filename,'time',time);
ncwrite(filename,'lat',lat_T66);
ncwrite(filename,'lon',lon_T66);
ncwrite(filename,'zi',z_T66_ref);
ncwrite(filename,'T_ice',round(T_fy_int_sur,2));
ncwrite(filename,'S_ice',round(S_fy_int_sur,2));
ncwrite(filename,'vb',round(vb_fy_int,0)/10);
ncwrite(filename,'rho',round(rho_fy_int_sur,1));
ncwrite(filename,'hi',-zT_fy_bot_sur);

ncwriteatt(filename,"/","title","MOSAiC first-year ice coring");
ncwriteatt(filename,"/","Conventions","CF-1.7");
ncwriteatt(filename,"/","standard_name_vocabulary","___");
ncwriteatt(filename,"/","contributor_name","Evgenii Salganik");
ncwriteatt(filename,"/","contributor_email","evgenii.salganik@proton.me");
ncwriteatt(filename,"/","institution","Norwegian Polar Institute (NPI)");
ncwriteatt(filename,"/","creator_name","Evgenii Salganik");
ncwriteatt(filename,"/","creator_email","evgenii.salganik@proton.me");
ncwriteatt(filename,"/","project","MOSAiC Single Column Modeling Working Group");
ncwriteatt(filename,"/","summary","First-year ice temperature, salinity, brine volume, and density");
ncwriteatt(filename,"/","id","DOI: XXXXXXXX");
ncwriteatt(filename,"/","license","CC-0");
ncwriteatt(filename,"/","metadata_link","http://www.doi-metatadat-link.org/mydoiiiii");
ncwriteatt(filename,"/","references","First-year sea-ice salinity, temperature, density, nutrient, oxygen and hydrogen isotope composition from the main coring site (MCS-FYI) during MOSAiC legs 1 to 4 in 2019/2020, version 2 https://doi.org/10.1594/PANGAEA.971385, Temperature and heating induced temperature difference measurements from SIMBA-type sea ice mass balance buoy 2019T66, deployed during MOSAiC 2019/20, https://doi.org/10.1594/PANGAEA.938134");
ncwriteatt(filename,"/","time_coverage_start","2019-10-29 06:00:00");
ncwriteatt(filename,"/","time_coverage_end","2020-07-29 06:30:00");
ncwriteatt(filename,"/","naming_authority","___");
ncwriteatt(filename,"/","keywords","surface energy budget, arctic, polar, salinity, temperature");
ncwriteatt(filename,"/","calendar","standard");
ncwriteatt(filename,"/","date_created","2025-01-02 15:00:00");
ncwriteatt(filename,"/","featureType","timeseries");

%% SYI coring, refined temperatures
clear; clc; 
load('T62.mat',"t_T62","T_T62","lat_T62","long_T62"); % Thermistor 2019T62 data import (SYI)
for i = 899:904; T_T62(i,:) = T_T62(898,:); end % remove errors
for i = 927:940; T_T62(i,:) = T_T62(926,:); end % remove errors
for i = 943; T_T62(i,:) = T_T62(942,:); end % remove errors
for i = 947:971; T_T62(i,:) = T_T62(946,:); end % remove errors
cut = 1062; t_T62 = t_T62(1:cut); T_T62 = T_T62(1:cut,:); lat_T62 = lat_T62(1:cut); lon_T62 = long_T62(1:cut);
td_T62 = datenum(t_T62); td_0_T62 = (datenum(t_T62)- datenum(t_T62(1))); % time of experiment [d]
top_t_dir_T62 = [0 212 225 234 254 272]; % Direct time of measurements [d]
top_dir_T62 =  -[0   0   0   0  36  52]/100; % Thickness for interpolation [m]
bot_t_dir_T62 = [0   5  46  81 129 185 220 230 266]; % Direct time of measurements [d]
bot_dir_T62 =  -[98 98 106 126 152 178 178 174 170]/100; % Thickness for interpolation [m]
bot_int_T62 = interp1(bot_t_dir_T62,bot_dir_T62,td_0_T62,'pchip'); % Experimental bottom interpolation
top_int_T62 = interp1(top_t_dir_T62,top_dir_T62,td_0_T62,'pchip'); % Experimental bottom interpolation
z_T62 = 0:-0.02:-4.78; z_T62 = z_T62 +0.92; % 2019T62 freeboard
z_T62_ref = 0.0:-0.005:-1.8; % IMB vertical coordinates, refined, 5 mm spacing
for i = 1:length(T_T62)
    T_T62_ref(i,:) = interp1(z_T62,T_T62(i,:),z_T62_ref,'linear'); % refinement of IMB temperature
end
frz_T62_ref = T_T62_ref; % temperature of only frozen sesnors
for i = 1:size(T_T62_ref,1) % time
    for j = 1:size(T_T62_ref,2) % depth
        if  j > 1+(z_T62_ref(1) - bot_int_T62(i))/-diff(z_T62_ref(1:2))
            frz_T62_ref(i,j) = NaN;
        elseif  j < 1+(z_T62_ref(1) - top_int_T62(i))/-diff(z_T62_ref(1:2))
            frz_T62_ref(i,j) = NaN;
        end
    end
end
clearvars err td_0_T62 top_t_dir_T62 top_dir_T62 bot_t_dir_T62 bot_dir_T62 T_T62 i j long_T62

% figure
% range = [-30 -25 -20 -15 -10 -8 -6 -4 -2 0 2 4]; 
% contourf(datenum(t_T62),z_T62_ref,frz_T62_ref',range,'-','ShowText','on','LabelSpacing',400,'LineColor','none','LineWidth',0.1); hold on
% plot(datenum(t_T62),top_int_T62,'k--'); plot(datenum(t_T62),bot_int_T62,'k--');
% hYLabel = ylabel('Depth (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
% clim([-30 5]);
% yticks(-2.0:0.2:0.2); ylim([-2.0 0.2]); % depth limits
% colorbar; brighten(.1); colormap(gca,bluewhitered);
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% hBar1 = colorbar; ylabel(hBar1,'Ice temperature (°C)','FontSize',8);
% clearvars range hYLabel ax hBar1

%% Stretch of T from IMB towards surface interface
T_sy_int = frz_T62_ref'; zT_sy_int = z_T62_ref;
for i = 1:length(T_sy_int)
     zT_sy_top(i) = zT_sy_int(find(~isnan(T_sy_int(:,i)),1));
     zT_sy_bot(i) = zT_sy_int(find(~isnan(T_sy_int(:,i)),1,'last'));
     zT_sy_top_sur(i) = zT_sy_top(i)-zT_sy_top(i);
     zT_sy_bot_sur(i) = zT_sy_bot(i)-zT_sy_top(i);
end
zT_sy_str = zT_sy_int(:)*ones(1,length(T_sy_int));
for i = 1:length(T_sy_int)
    zT_sy_str(:,i) = (zT_sy_str(:,i)-zT_sy_top(i))*(+zT_sy_bot_sur(i)-zT_sy_top_sur(i))/(zT_sy_bot(i)-zT_sy_top(i)) + zT_sy_top_sur(i);
end
T_sy_int_sur = T_sy_int; % preallocation
for i = 1:length(T_sy_int)
    T_sy_int_sur(:,i) = interp1(zT_sy_str(:,i),T_sy_int(:,i),zT_sy_int','linear'); % interpolation of IMB temp along salinity core depth
end
[X,Y] = meshgrid(datenum(t_T62),zT_sy_int); Z = T_sy_int_sur; % contour plot preparation
clearvars i

figure
range = [-30 -25 -20 -15 -10 -8 -6 -4 -2 0 2 4]; 
contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','none','LineWidth',0.1); hold on
plot(datenum(t_T62),zT_sy_top-zT_sy_top,'k--'); plot(datenum(t_T62),zT_sy_bot-zT_sy_top,'k--');
hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
clim([-30 5]);
yticks(-2.0:0.2:0.2); ylim([-2.0 0.2]); % depth limits
colorbar; brighten(.1); colormap(gca,bluewhitered);
ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
hBar1 = colorbar; ylabel(hBar1,'Ice temperature (°C)','FontSize',8);
clearvars X Y Z range hYLabel ax hBar1

%% Salinity base data
load('Coring_old_fb.mat',"S_sy","zS_sy","zzS_sy","t_sy","td_sy","td"); % Coring data import
for i = 1:length(t_sy)
    [~,idx(i)] = min(abs(datenum(t_sy(i))-datenum(t_T62)));
    l_imb(i) = -zT_sy_bot(idx(i))+zT_sy_top(idx(i)); % core length from IMB for coring dates
end
dS_sy = zS_sy;
for i = 1:length(S_sy)
    hS(i) = zzS_sy{i,1}(end,2); % salinity core length
    dS_sy{i,1} = -zS_sy{i,1}*l_imb(i)/hS(i); % length of salinity measurements for FYI
end
S_lvl = S_sy; S_lvl(17)=[]; S_lvl(13)=[]; S_lvl(12)=[];
dS_sy_lvl = dS_sy; dS_sy_lvl(17)=[]; dS_sy_lvl(13)=[]; dS_sy_lvl(12)=[];
td_lvl = td_sy; td_lvl(17)=[]; td_lvl(13)=[]; td_lvl(12)=[];
S_vect = cell2mat(S_lvl); dS_sy_vect = cell2mat(dS_sy_lvl); tdS_sy_vect = cell2mat(td_lvl);
x = tdS_sy_vect; y = dS_sy_vect; z = S_vect;

xv = datenum(t_T62);
yv = z_T62_ref;
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y); S_sy_int = Z;
clearvars l_imb lvl tdS_sy_vect dS_sy_vect S_vect idx x y z xv yv i S_lvl td td_sy t_sy hS td_lvl dS_sy_lvl zzS_sy

% figure
% range = [0 0.2 1 3 5 6 8 10]; % Temperature graph accuracy;
% [C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
% plot(datenum(t_T62),zT_sy_top-zT_sy_top,'k--'); plot(datenum(t_T62),zT_sy_bot-zT_sy_top,'k--');
% clabel(C,h,'FontSize',6,'Color','k');
% hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); colorbar; colormap summer; brighten(.7);
% clim([0 8]);
% yticks(-2.0:0.2:0.2); ylim([-2.0 0.2]);
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% hBar1 = colorbar; ylabel(hBar1,'Ice salinity','FontSize',8);
% clearvars X Y Z C h range hYLabel ax hBar1

%% Salinity stretch over IMB depth
for i = 1:length(S_sy_int)
     zS_sy_top(i) = z_T62_ref(find(~isnan(S_sy_int(:,i)),1));
     zS_sy_bot(i) = z_T62_ref(find(~isnan(S_sy_int(:,i)),1,'last'));
     zT_sy_top(i) = z_T62_ref(find(~isnan(T_sy_int_sur(:,i)),1));
     zT_sy_bot(i) = z_T62_ref(find(~isnan(T_sy_int_sur(:,i)),1,'last'));
end
zT_sy_str = z_T62_ref(:)*ones(1,length(T_sy_int));
for i = 1:length(T_sy_int)
    zT_sy_str(:,i) = (zT_sy_str(:,i)-zS_sy_top(i))*(+zT_sy_bot(i)-zT_sy_top(i))/(zS_sy_bot(i)-zS_sy_top(i)) + zT_sy_top(i);
end

S_sy_int_sur = T_sy_int; % preallocation
for i = 1:length(T_sy_int)
    S_sy_int_sur(:,i) = interp1(zT_sy_str(:,i),S_sy_int(:,i),z_T62_ref','linear'); % interpolation of salinity along IMB depth
end
[X,Y] = meshgrid(datenum(t_T62),z_T62_ref); Z = S_sy_int_sur; % contour plot preparation
clearvars zT_sy_str i

figure
range = [0 0.2 1 2 3 4 6 8 10]; % Salinity graph accuracy;
[C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
plot(datenum(t_T62),zT_sy_top-zT_sy_top,'k--'); plot(datenum(t_T62),zT_sy_bot-zT_sy_top,'k--');
clabel(C,h,'FontSize',7,'Color','k');
hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
clim([0 8]);
yticks(-2.0:0.2:0.2); ylim([-2.0 0.2]);
colorbar; colormap summer; brighten(.7);
ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
hBar1 = colorbar; ylabel(hBar1,'Ice salinity','FontSize',8);
clearvars X Y Z C range h hYLabel ax hBar1

%% Brine volume
vb_sy_int = zeros(size(S_sy_int_sur)); F1_sy_int = vb_sy_int; F2_sy_int = vb_sy_int; rhoi_sy_int = vb_sy_int;
for j = 1:length(t_T62)
    for i = 1:length(z_T62_ref)
        if T_sy_int_sur(i,j) <= -2
            F1_sy_int(i,j) = -4.732-22.45*T_sy_int_sur(i,j) - 0.6397*T_sy_int_sur(i,j).^2 - 0.01074*T_sy_int_sur(i,j).^3; % Cox and Weeks (1983)
            F2_sy_int(i,j) = 8.903*10^-2 - 1.763*10^-2*T_sy_int_sur(i,j) - 5.33*10^-4*T_sy_int_sur(i,j)^2 - 8.801*10^-6*T_sy_int_sur(i,j)^3; % Cox and Weeks (1983)
        else
            F1_sy_int(i,j) = -4.1221*10^-2 + -1.8407*10^1*T_sy_int_sur(i,j) + 5.8402*10^-1*T_sy_int_sur(i,j)^2 + 2.1454*10^-1*T_sy_int_sur(i,j)^3; % F1 from Lepparanta and Manninen (1988)
            F2_sy_int(i,j) = 9.0312*10^-2 + -1.6111*10^-2*T_sy_int_sur(i,j) + 1.2291*10^-4*T_sy_int_sur(i,j)^2 + 1.3603*10^-4*T_sy_int_sur(i,j)^3; % F2 from Lepparanta and Manninen (1988)
        end
        rhoi_sy_int(i,j) = (917-1.403*10^-1*T_sy_int_sur(i,j)); % pure ice density, Pounder (1965) for T_insitu
        vb_sy_int(i,j) = (1000).*(rhoi_sy_int(i,j).*S_sy_int_sur(i,j)/1000)./(F1_sy_int(i,j)-rhoi_sy_int(i,j).*S_sy_int_sur(i,j)/1000.*F2_sy_int(i,j)); % Brine volume for T_insitu, LM
        if vb_sy_int(i,j) < 0
            vb_sy_int(i,j) = NaN;
        elseif vb_sy_int(i,j) > 1000
            vb_sy_int(i,j) = NaN;
        end
    end
end
clearvars i j F1_sy_int F2_sy_int rhoi_sy_int

% figure
% [X,Y] = meshgrid(datenum(t_T62),z_T62_ref); Z = vb_sy_int/1000; % contour plot preparation
% range = [0 0.02 0.05 0.1 0.2 0.3];
% [C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w','LineWidth',0.1); hold on
% plot(datenum(t_T62),zT_sy_top-zT_sy_top,'k--'); plot(datenum(t_T62),zT_sy_bot-zT_sy_top,'k--');
% clabel(C,h,'FontSize',7,'Color','k');
% hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
% colorbar; colormap summer; brighten(0.5); clim([0 0.3]);
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% yticks(-2.2:0.2:0.4); ylim([-2.0 0.2]); % summertime depth limits
% hBar1 = colorbar; ylabel(hBar1,'Brine volume','FontSize',8);
% clearvars X Y Z range hYLabel ax hBar1

%% Density base data
load('Coring_old_fb.mat',"rho_sy","zrho_sy","td_rho_sy","t_sy"); % Coring data import
for i = [1:3 5:18]
    [~,idx(i)] = min(abs(datenum(t_sy(i))-datenum(t_T62)));
    l_imb(i) = -zT_sy_bot(idx(i))+zT_sy_top(idx(i)); % core length from IMB for coring dates
end
for i = [1:3 5:18]
    h_rho(i) = zrho_sy{i,1}(end); % density core length
    zrho_sy{i,1} = -zrho_sy{i,1}*l_imb(i)/h_rho(i); % length of density measurements for SYI
end
rho_lvl = rho_sy; drho_sy_lvl = zrho_sy; td_lvl = td_rho_sy;
drho_sy_lvl{2}(2)=[]; rho_lvl{2}(2)=[]; td_lvl{2}(2)=[]; % remove NaN week 2
rho_vect = cell2mat(rho_lvl); drho_sy_vect = cell2mat(drho_sy_lvl); tdvg_sy_vect = cell2mat(td_lvl);
x = tdvg_sy_vect; y = drho_sy_vect; z = rho_vect;
xv = datenum(t_T62); yv = z_T62_ref;
[X,Y] = meshgrid(xv,yv);
Z = griddata(x,y,z,X,Y); rho_sy_int = Z;
clearvars l_imb lvl tdS_sy_vect dS_sy_vect S_vect idx x y z xv yv i S_lvl td td_sy t_sy hS td_lvl dS_sy_lvl zzS_sy

% figure
% range = [700 750 800 840 860 900 910 940]; % density graph accuracy
% [C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
% plot(datenum(t_T62),zT_sy_top-zT_sy_top,'k--'); plot(datenum(t_T62),zT_sy_bot-zT_sy_top,'k--');
% clabel(C,h,'FontSize',7,'Color','k');
% hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
% colorbar; colormap summer; brighten(0.5);
% clim([800 950]);
% yticks(-2:0.2:0.2); ylim([-2 0.2]); % depth limits
% ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
% hBar1 = colorbar; ylabel(hBar1,'Ice density (kg/m^3)','FontSize',8);

%% Density stretch over IMB depth
for i = 1:length(rho_sy_int)
     zS_sy_top(i) = z_T62_ref(find(~isnan(rho_sy_int(:,i)),1));
     zS_sy_bot(i) = z_T62_ref(find(~isnan(rho_sy_int(:,i)),1,'last'));
     zT_sy_top(i) = z_T62_ref(find(~isnan(T_sy_int_sur(:,i)),1));
     zT_sy_bot(i) = z_T62_ref(find(~isnan(T_sy_int_sur(:,i)),1,'last'));
end
zT_sy_str = z_T62_ref(:)*ones(1,length(T_sy_int));
for i = 1:length(T_sy_int)
    zT_sy_str(:,i) = (zT_sy_str(:,i)-zS_sy_top(i))*(+zT_sy_bot(i)-zT_sy_top(i))/(zS_sy_bot(i)-zS_sy_top(i)) + zT_sy_top(i);
end

rho_sy_int_sur = T_sy_int; % preallocation
for i = 1:length(T_sy_int)
    rho_sy_int_sur(:,i) = interp1(zT_sy_str(:,i),rho_sy_int(:,i),z_T62_ref','linear'); % interpolation of salinity along IMB depth
end
[X,Y] = meshgrid(datenum(t_T62),z_T62_ref); Z = rho_sy_int_sur; % contour plot preparation
clearvars zT_sy_str i

figure
range = [700 750 800 840 860 900 910 940]; % density graph accuracy
[C,h] = contourf(X,Y,Z,range,'-','ShowText','on','LabelSpacing',400,'LineColor','w'); hold on
plot(datenum(t_T62),zT_sy_top-zT_sy_top,'k--'); plot(datenum(t_T62),zT_sy_bot-zT_sy_top,'k--');
clabel(C,h,'FontSize',6,'Color','k');
hYLabel = ylabel('Ice thickness (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
clim([800 950]);
yticks(-2.0:0.2:0.2); ylim([-2.0 0.2]); % depth limits
colorbar; colormap summer; brighten(0.5);
ax = gca; ax.XTick = datenum(datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
hBar1 = colorbar; ylabel(hBar1,'Ice density (kg/m^3)','FontSize',8);
clearvars X Y Z C range h hYLabel ax hBar1

%% netCDF export
filename = 'MOSAiC_SYI_coring.nc';
time = datenum(t_T62-datetime('1979-01-01 00:00:00')); % time to nd datenum
delete(filename)

nccreate(filename,'time','Dimensions',{'time' length(time)},'format','netcdf4'); % time
ncwriteatt(filename,'time','standard_name','time');
ncwriteatt(filename,'time','long_name','time');
ncwriteatt(filename,'time','units','days since 1979-01-01 00:00:00');

nccreate(filename,'lat','Dimensions',{'time' length(time)});
ncwriteatt(filename,'lat','standard_name','lat');
ncwriteatt(filename,'lat','long_name','latitude');
ncwriteatt(filename,'lat','units','°');

nccreate(filename,'lon','Dimensions',{'time' length(time)});
ncwriteatt(filename,'lon','standard_name','lon');
ncwriteatt(filename,'lon','long_name','longitude');
ncwriteatt(filename,'lon','units','°');

nccreate(filename,'zi','Dimensions',{'zi' length(z_T62_ref)}); % Ice core length
ncwriteatt(filename,'zi','standard_name','zi');
ncwriteatt(filename,'zi','long_name','Ice core length');
ncwriteatt(filename,'zi','units','m');

nccreate(filename,'T_ice','Dimensions',{"zi",length(z_T62_ref),"time",length(time)}); % Temperature
ncwriteatt(filename,'T_ice','standard_name','Temperature');
ncwriteatt(filename,'T_ice','long_name','Ice temperature from 2019T62 buoy');
ncwriteatt(filename,'T_ice','units','°C');

nccreate(filename,'S_ice','Dimensions',{"zi",length(z_T62_ref),"time",length(time)}); % Salinity
ncwriteatt(filename,'S_ice','standard_name','Salinity');
ncwriteatt(filename,'S_ice','long_name','Ice practical salinity');
ncwriteatt(filename,'S_ice','units','PSU');

nccreate(filename,'vb','Dimensions',{"zi",length(z_T62_ref),"time",length(time)}); % Brine volume
ncwriteatt(filename,'vb','standard_name','Brine volume');
ncwriteatt(filename,'vb','long_name','Brine volume fraction');
ncwriteatt(filename,'vb','units','%');

nccreate(filename,'rho','Dimensions',{"zi",length(z_T62_ref),"time",length(time)}); % Ice density
ncwriteatt(filename,'rho','standard_name','Ice density');
ncwriteatt(filename,'rho','long_name','Ice density at -15°C laboratory temperature');
ncwriteatt(filename,'rho','units','kg/m3');

nccreate(filename,'hi','Dimensions',{'time' length(time)}); % Ice thickness
ncwriteatt(filename,'hi','standard_name','hi');
ncwriteatt(filename,'hi','long_name','Ice thickness from 2019T62 buoy');
ncwriteatt(filename,'hi','units','m');

ncwrite(filename,'time',time);
ncwrite(filename,'lat',lat_T62);
ncwrite(filename,'lon',lon_T62);
ncwrite(filename,'zi',z_T62_ref);
ncwrite(filename,'T_ice',round(T_sy_int_sur,2));
ncwrite(filename,'S_ice',round(S_sy_int_sur,2));
ncwrite(filename,'vb',round(vb_sy_int,0)/10);
ncwrite(filename,'rho',round(rho_sy_int_sur,1));
ncwrite(filename,'hi',-zT_sy_bot_sur);

ncwriteatt(filename,"/","title","MOSAiC second-year ice coring");
ncwriteatt(filename,"/","Conventions","CF-1.7");
ncwriteatt(filename,"/","standard_name_vocabulary","___");
ncwriteatt(filename,"/","contributor_name","Evgenii Salganik");
ncwriteatt(filename,"/","contributor_email","evgenii.salganik@proton.me");
ncwriteatt(filename,"/","institution","Norwegian Polar Institute (NPI)");
ncwriteatt(filename,"/","creator_name","Evgenii Salganik");
ncwriteatt(filename,"/","creator_email","evgenii.salganik@proton.me");
ncwriteatt(filename,"/","project","MOSAiC Single Column Modeling Working Group");
ncwriteatt(filename,"/","summary","Second-year ice temperature, salinity, brine volume, and density");
ncwriteatt(filename,"/","id","DOI: XXXXXXXX");
ncwriteatt(filename,"/","license","CC-0");
ncwriteatt(filename,"/","metadata_link","http://www.doi-metatadat-link.org/mydoiiiii");
ncwriteatt(filename,"/","references","Second-year sea-ice salinity, temperature, density, nutrient, oxygen and hydrogen isotope composition from the main coring site (MCS-SYI) during MOSAiC legs 1 to 4 in 2019/2020, https://doi.pangaea.de/10.1594/PANGAEA.973860, Temperature and heating induced temperature difference measurements from SIMBA-type sea ice mass balance buoy 2019T62, deployed during MOSAiC 2019/20, https://doi.org/10.1594/PANGAEA.940231");
ncwriteatt(filename,"/","time_coverage_start","2019-10-29 02:30:00");
ncwriteatt(filename,"/","time_coverage_end","2020-07-20 08:30:00");
ncwriteatt(filename,"/","naming_authority","___");
ncwriteatt(filename,"/","keywords","surface energy budget, arctic, polar, salinity, temperature");
ncwriteatt(filename,"/","calendar","standard");
ncwriteatt(filename,"/","date_created","2025-01-04 23:00:00");
ncwriteatt(filename,"/","featureType","timeseries");

%% import
close all; clear; clc; 
project = 'MOSAiC_FYI_coring.nc'; ncdisp(project);
time = ncread(project,'time'); t_0 = (datetime('1979-01-01 00:00:00')); t = t_0 + days(time);
zi = ncread(project,'zi'); T_ice = ncread(project,'T_ice'); S_ice = ncread(project,'S_ice'); vb = ncread(project,'vb'); rho = ncread(project,'rho'); hi = ncread(project,'hi');