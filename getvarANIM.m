%%%% This script will create animations of a chosen tracer for model outputs 
% Written: Jonathan Rogerson
% Date: June 2024
% Adapted for the FTLE computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNIVERSITY of CAPE TOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script calls CROCO outputs 
% Objective is to animate over variable of choice fo the model period
% Remember, visual outputs are very expensive 


addpath /home/jono/CROCO/croco_tools/UTILITIES/m_map1.4h
addpath /media/data/DATASETS_CROCOTOOLS/m_map1.4f
addpath /home/jono/Documents/MATLAB/CMOCEAN
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt'
addpath '/home/jono/Documents/MATLAB/CDT-master/cdt/cdt_data'
addpath '/home/jono/Documents/MATLAB/tight_subplot'
addpath '/home/jono/Documents/MATLAB/LCS-Tool-master'
addpath '/usr/local/MATLAB/R2020a/toolbox/matlab/imagesci/'
addpath '/home/jono/CROCO/croco_tools/Preprocessing_tools'

%% Get the main File paths for the input data

load FTLE_zoom_nonorm_2014.mat
time = length(mydate);
var = myFTLE;
%% Figure and animation

cmin = 0;
cmax = 0.25;
levels = [cmin:0.05:cmax];
%mydepths = [100,200,300];

for i = 259:time
    clf;
    clc;
    figure(1); 
    hold on;
    title(string(mydate(:,i)),'fontsize',15);

    % Plot loop
    m_proj('miller','long',[min(min(lon)), max(max(lon))], ...
                'lat',[min(min(lat)), max(max(lat))]);
    m_pcolor(lon,lat,var(:,:,i));
    shading flat
    cmocean('ice');
    hold on;
%    m_contour(CROCO_lon,CROCO_lat,CROCO_top,mydepths,'Color',...
%        [0.6824,0.9686,0.6118],'LineWidth',1,'LineStyle','--');
    m_gshhs_f('patch',[.8 .8 .8],'edgecolor','none');
    m_grid('box','fancy','linest','none','tickdir','out','fontsize',13);
    caxis([cmin cmax])
    cRange=caxis;
    ca = colorbar;
    ca.Label.String = 'FTLE [days^{-1}]';
    ca.FontSize = 12;
    caxis([cmin cmax]);
    
    % Animate and save
    pause(0.02);
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    filename = strcat('FTLE_SBUS_anim',string(i),'.gif');
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    
    %if i == 1
    %     imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    %else
    %     imwrite(imind,cm,filename,'gif','WriteMode','append');
    %end 
end








