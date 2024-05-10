%%%% This script will create a composite figure with the aim of showing the
%%%% tranpport stream function for mean summer and winter and also overlay
%%%% the FTLE field

% We already have some problems a the amount of data to process is huge. So
% to save space we will on the offset only read in the summer and winter
% files and process them appropriately. This will be done to save time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD all the neccesarily library paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/jono/CROCO/CONFIGS/SBUS_3km
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools/Transport
addpath /home/jono/CROCO/croco_tools/Diagnostic_tools/EKE
addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath /home/jono/Documents/MATLAB/PIFF-master
addpath '/home/jono/Documents/MATLAB/LCS-Tool-master'

% Run the start and crocotools to get varibles into the worksapce
start
%crocotools_param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
CROCO_path = '/media/data/CHPC_SBUS_3km';
Ymin = 2005;
Ymax = 2005;
Yorig = 1990;

%  !!! WARNING weak point: vtransform should be the one used for CROCO
%
vtransform=1;
%
% Lower level (NaN=bottom)
%
z1=-100;
%
% Upper level (NaN=surface)
%
z2=NaN;
%
% Output matlab file
%
outname='transport_croco.mat';

%%%%%%%%%
% CROCO
%%%%%%%%%

% Once-off gridding data and dimensions
file = strcat(CROCO_path,'/','croco_avg_Y',string(Ymin),'M',string(1),'.nc');
file = convertStringsToChars(file);

nc=netcdf(file);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
rmask=nc{'mask_rho'}(:);
h=nc{'h'}(:);
theta_s=nc.theta_s(:);
theta_b=nc.theta_b(:);
hc=nc.hc(:);
N=length(nc('s_rho'));
mask=rmask;
mask(mask==0)=NaN;
close(nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STREAM function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% SUMMER
%%%%%%%%%%%%%%%%%%%%%%%
summer = [1, 2, 12];  % Dec, Jan, Feb  

PSI_sum = [];
U_sum = [];
V_sum = [];
TIME_sum = [];

for i = Ymin:Ymax
    for j = summer
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        file = convertStringsToChars(file);
        disp(file)
        nc=netcdf(file);
        time = nc{'time'}(:);
        for k = 1:length(time)
            zeta=squeeze(nc{'zeta'}(k,:,:));
            u=squeeze(nc{'u'}(k,:,:,:));
            v=squeeze(nc{'v'}(k,:,:,:));
            %
        %  Get the transport between the surface and z0
            zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);
            zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
            %
            [u,hu]=vintegr2(u,rho2u_3d(zw),rho2u_3d(zr),z1,z2);
            [v,hv]=vintegr2(v,rho2v_3d(zw),rho2v_3d(zr),z1,z2);
            %
            % Compute PSI
            %
            [u,v]=get_obcvolcons(u,v,pm,pn,rmask,[1 1 1 1]);
            [psi0,psi1,island]=get_psi0(u,v,pm,pn,rmask);  
            if sum(sum(island))==0
                A=0;
            else
                A=get_a(u,v,psi0,psi1,island,pm,pn);
            end
            psi=psi0+A*psi1;
            % Store arrays
            PSI_sum = cat(3,PSI_sum,psi);
            U_sum = cat(3,U_sum,u);
            V_sum = cat(3,V_sum,v);
            % Memory issue
            clear u v zeta zw zr hu hv psi0 psi1 island
        end
        TIME_sum = cat(1,TIME_sum,time);
        close(nc);
        disp('Stored Month PSI and time')
    else
        disp(strcat('No data for',file))
    end
    end
end

clear time u v 

% Clean-up

[~, ia, ~] = unique(TIME_sum);
PSI_sum = PSI_sum(:,:,ia);
U_sum = U_sum(:,:,ia);
V_sum = V_sum(:,:,ia);
TIME_sum = TIME_sum(ia);

% Save

%save(strcat('summer_',outname),'lon','lat','mask','TIME_sum','U_sum',...
%    'V_sum','PSI_sum','h','z1','z2')

%%%%%%%%%%%%%%%%%%%%%%%
% WINTER
%%%%%%%%%%%%%%%%%%%%%%%

winter = [6, 7, 8]; 

PSI_win = [];
U_win = [];
V_win = [];
TIME_win = [];

for i = Ymin:Ymax
    for j = winter
    file = strcat(CROCO_path,'/','croco_avg_Y',string(i),'M',string(j),'.nc');
    if exist(file, 'file')
        file = convertStringsToChars(file);
        disp(file)
        nc=netcdf(file);
        time = nc{'time'}(:);
        for k = 1:length(time)
            zeta=squeeze(nc{'zeta'}(k,:,:));
            u=squeeze(nc{'u'}(k,:,:,:));
            v=squeeze(nc{'v'}(k,:,:,:));
            %
        %  Get the transport between the surface and z0
            zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);
            zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
            %
            [u,hu]=vintegr2(u,rho2u_3d(zw),rho2u_3d(zr),z1,z2);
            [v,hv]=vintegr2(v,rho2v_3d(zw),rho2v_3d(zr),z1,z2);
            %
            % Compute PSI
            %
            [u,v]=get_obcvolcons(u,v,pm,pn,rmask,[1 1 1 1]);
            [psi0,psi1,island]=get_psi0(u,v,pm,pn,rmask);  
            if sum(sum(island))==0
                A=0;
            else
                A=get_a(u,v,psi0,psi1,island,pm,pn);
            end
            psi=psi0+A*psi1;
            % Store arrays
            PSI_win = cat(3,PSI_win,psi);
            U_win = cat(3,U_win,u);
            V_win = cat(3,V_win,v);
            % Memory issue
            clear u v zeta zw zr hu hv psi0 psi1 island
        end
        TIME_win = cat(1,TIME_win,time);
        close(nc);
        disp('Stored Month PSI and time')
    else
        disp(strcat('No data for',file))
    end
    end
end

clear time u v 

% Clean-up

[~, ia, ~] = unique(TIME_win);
PSI_win = PSI_win(:,:,ia);
U_win = U_win(:,:,ia);
V_win = V_win(:,:,ia);
TIME_win = TIME_win(ia);

%save(strcat('winter_',outname),'lon','lat','mask','TIME_win','U_win',...
%    'V_win','PSI_win','h','z1','z2')

%% Divide between summer and winter
% Convert the time data

TIME = sort([TIME_sum; TIME_win]);

CROCO_time = datetime(Yorig,1,1) + seconds(TIME);
[Y,MO,D] = datevec(CROCO_time);  

% Average all variable data for our seasons

% Index U, V and psi

U_sum = double(nanmean(U_sum,3));
V_sum = double(nanmean(V_sum,3));

U_win = double(nanmean(U_win,3));
V_win = double(nanmean(V_win,3));

PSI_sum = double(nanmean(PSI_sum,3));
PSI_win = double(nanmean(PSI_win,3));

%% Plot the figure

figure
npts=[2 2 2 2];

lonmin=15;
lonmax=20;
latmin=-36;
latmax=-28;

topo=1;

subplot(1,2,1)  % Summer

m_proj('mercator',...
       'lon',[lonmin lonmax],...
       'lat',[latmin latmax]);
[x,y]=m_ll2xy(lon,lat,'clip','off');
hold on

psi_r=1e-6*psi2rho(PSI_sum).*mask;
u = U_sum;
v = V_sum;

colmax = 46;
colmin = -12;
dcol = 8;

[C1,h1]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[dcol:dcol:colmax],'k');
if ~isempty(C1)
  clabel(C1,h1,'LabelSpacing',1000,'Rotation',0)
  hf1=add_streamarrows(C1,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C2,h2]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(-psi_r,npts)...
           ,[dcol:dcol:colmin],'k');
if ~isempty(C2)
  clabel(C2,h2,'LabelSpacing',1000,'Rotation',0)
  hf2=add_streamarrows(C2,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C3,h3]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[0 0],'k');
if ~isempty(C3)
  clabel(C3,h3,'LabelSpacing',1000,'Rotation',0)
  hf3=add_streamarrows(C3,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h3,'Color','k','Linewidth',1.2);
  set(hf3,'Linewidth',1.2);
  hold on
end

[C4,h4]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[22:1:26],'k');
if ~isempty(C4)
  clabel(C4,h4,'LabelSpacing',1000,'Rotation',0)
  hf4=add_streamarrows(C4,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h4,'Color','k','Linewidth',1.2);
  set(hf4,'Linewidth',1.2);
  hold on
end

if topo==1
  [C5,h5]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(h,npts)...
           ,[100 200 500],'--g');
  %clabel(C5,h5,'LabelSpacing',1000,'Rotation',0)
end

m_usercoast('coastline_l.mat','patch',[.9 .9 .9]);
hold off
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
set(findobj('tag','m_grid_color'),'facecolor','white')

if isfinite(z1)
 zname1=[num2str(abs(z1)),' m'];
else
 zname1='fond';
end
if isfinite(z2)
 zname2=[num2str(abs(z2)),' m'];
else
 zname2='surface';
end
title(['Transport summer [Svd] ',zname2,' - ',zname1])
hold off

subplot(1,2,2)

m_proj('mercator',...
       'lon',[lonmin lonmax],...
       'lat',[latmin latmax]);
[x,y]=m_ll2xy(lon,lat,'clip','off');
hold on

psi_r=1e-6*psi2rho(PSI_win).*mask;
u = U_win;
v = V_win;

colmax = round(max(max(psi_r)))+1;  % 11 winter
colmin = round(min(min(psi_r)))-1;  % -7 winter
dcol=1;
[C1,h1]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[dcol:dcol:colmax],'k');
if ~isempty(C1)
  clabel(C1,h1,'LabelSpacing',1000,'Rotation',0)
  hf1=add_streamarrows(C1,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C2,h2]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(-psi_r,npts)...
           ,[dcol:dcol:colmin],'k');
if ~isempty(C2)
  clabel(C2,h2,'LabelSpacing',1000,'Rotation',0)
  hf2=add_streamarrows(C2,x,y,u2rho_2d(u),v2rho_2d(v));
  hold on
end

[C3,h3]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[0 0],'k');
if ~isempty(C3)
  clabel(C3,h3,'LabelSpacing',1000,'Rotation',0)
  hf3=add_streamarrows(C3,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h3,'Color','k','Linewidth',1.2);
  set(hf3,'Linewidth',1.2);
  hold on
end

[C4,h4]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(psi_r,npts)...
           ,[20:0.1:21],'k');
if ~isempty(C4)
  clabel(C4,h4,'LabelSpacing',1000,'Rotation',0)
  hf4=add_streamarrows(C4,x,y,u2rho_2d(u),v2rho_2d(v));
  set(h4,'Color','k','Linewidth',1.2);
  set(hf4,'Linewidth',1.2);
  hold on
end

if topo==1
  [C5,h5]=m_contour(rempoints(lon,npts),rempoints(lat,npts),rempoints(h,npts)...
           ,[100 200 500],'--g');
  %clabel(C5,h5,'LabelSpacing',1000,'Rotation',0)
end

m_usercoast('coastline_l.mat','patch',[.9 .9 .9]);
hold off
m_gshhs_h('patch',[.8 .8 .8],'edgecolor','none');
m_grid('box','fancy','xtick',5,'ytick',5,'tickdir','out');
set(findobj('tag','m_grid_color'),'facecolor','white')

if isfinite(z1)
 zname1=[num2str(abs(z1)),' m'];
else
 zname1='fond';
end
if isfinite(z2)
 zname2=[num2str(abs(z2)),' m'];
else
 zname2='surface';
end
title(['Transport winter [Svd] ',zname2,' - ',zname1])
hold off