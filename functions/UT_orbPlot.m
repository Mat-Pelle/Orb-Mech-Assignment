% close all
% clear
% clc
% 
% addpath 'functions'
% addpath 'functions/time'

%% Parameters & initial conditions----------------------------------------
Earth_mu = UT_astroConstants(13);
Moon_mu  = UT_astroConstants(20);

j2      = UT_astroConstants(9);
Earth_R = UT_astroConstants(23);

orb_kep0 = [34814, 0.5054,deg2rad(51.6177),deg2rad(156.0701),deg2rad(253.1633),0];
orb_T = 2*pi*sqrt((orb_kep0(1)^3)/Earth_mu);

[orb_r0, orb_v0] = UT_kep2car(orb_kep0,Earth_mu);
orb_y0 = [orb_r0; orb_v0]';

date0 = [2027 4 1 0 0 0];
date0 = date2mjd2000(date0);

%% ---Options for the ODE solver & Timespan
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

orbPlot_n          = 16*500;                     %number of orbits to simulate
tmax               = orbPlot_n*orb_T;              %max time of simulation (nr of orbits*orbit period), in seconds.

orbPlot_stepsn     = 100;                      %number of segments per full orbit
orbPlot_stepst     = orb_T/orbPlot_stepsn;         %time per step
orbPlot_stepsnTot  = orbPlot_stepsn*orbPlot_n;         %total number of steps for the simulation

orbPlot_stepsYrs   = floor(tmax/(3600*24*365));
orbPlot_stepsDays  = floor(tmax/(3600*24));
orbPlot_stepsPerYr = floor((3600*24*365)/orbPlot_stepst);

orbPlot_tspan = linspace(0, tmax, orbPlot_stepsnTot);

%%
[orbPlot_t, orbPlot_car] = ode113(@(t,y)... 
    orbProp_carJ2Moon(t,y,Earth_mu,Moon_mu,j2,Earth_R,date0),...
    orbPlot_tspan, orb_y0, options);
%orbPlot_r = [orbPlot_car(:,1) orbPlot_car(:,2) orbPlot_car(:,3)];
%orbPlot_v = [orbPlot_car(:,4) orbPlot_car(:,5) orbPlot_car(:,6)];

%% Plotting

orbPlt_cMap = parula(orbPlot_stepsYrs);
orbPlt_year = 0;

figure(10)    
for j=1:orbPlot_stepsYrs
    
    orbPlt_color = orbPlt_cMap(j,:);
    
    orbPlt_i(j) = 1+orbPlot_stepsPerYr*orbPlt_year;
    plot3(orbPlot_car(orbPlt_i(j):orbPlt_i(j)+orbPlot_stepsn,1),...
          orbPlot_car(orbPlt_i(j):orbPlt_i(j)+orbPlot_stepsn,2),...
          orbPlot_car(orbPlt_i(j):orbPlt_i(j)+orbPlot_stepsn,3),...
          'HandleVisibility', 'off','LineWidth', 1.4,'Color',orbPlt_color);
    hold on
    orbPlt_year = orbPlt_year+1;
end

orbPlot_c = colorbar('Ticks',[-6000,0,6000],...
         'TickLabels',{'0','8','16'});
orbPlot_c.Label.String = 'Delta t [years]';

hold on
[Earth_x,Earth_y,Earth_z]=sphere(180);
Earth_x = Earth_x*Earth_R;
Earth_y = Earth_y*Earth_R;
Earth_z = Earth_z*Earth_R;
Earth = surf(Earth_x, Earth_y, Earth_z,'EdgeColor','none','FaceAlpha',0.1);
map=imread('EarthTexture.jpg');
map=imrotate(map,180);
warp(Earth_x,Earth_y,Earth_z,map);
axis equal
grid on
title('Orbit evolution over time');