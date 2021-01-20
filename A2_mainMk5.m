%% Startup----------------------------------------------------------------
close all
clear
clc

addpath 'functions'
addpath 'functions/time'

%% Parameters & initial conditions----------------------------------------

% constants --------------------------------------------------------------
Earth_mu = UT_astroConstants(13);
Moon_mu  = UT_astroConstants(20);

j2      = UT_astroConstants(9);
Earth_R = UT_astroConstants(23);

%Initial conditions ------------------------------------------------------
%Keplerian elements
orb_kep0 = [34814, 0.5054,deg2rad(51.6177),deg2rad(156.0701),deg2rad(253.1633),0];
orb_T = 2*pi*sqrt((orb_kep0(1)^3)/Earth_mu);

%Cartesian coordinates
[orb_r0, orb_v0] = UT_kep2car(orb_kep0,Earth_mu);
orb_y0 = [orb_r0; orb_v0]';

%Start date
date0 = [2027 4 1 0 0 0];
date0 = date2mjd2000(date0);

sim_n = 4; %'1'=groundtrack '2'=orbit propagation '3'=orbit evolution '4'=real data comparison
switch sim_n

case 1 % Groundtrack------------------------------------------------------
    
    %Options for the ODE solver & Timespan--------------------------------
    gt_we        = deg2rad(15.04/3600);   %Angular speed of Earth [rad/s]
    gt_longGreen = 0;                     %Greenwich longitude [rad]

    %Different cases
    n = 1;  %input: '1' = 1 orbit, '2' = 1 day, '3' = 10 days
    switch n
    case 1
        gt_tmax = orb_T;
    
    case 2
        gt_tmax = 86400;
    
    case 3
        gt_tmax = 864000;

    end 

    %Set time span
    gt_nOrb = gt_tmax/orb_T;
    gt_t0 = 0;
    gt_tspan = linspace(gt_t0,gt_tmax,100000);

    %Groundtrack calculations-------------------------------------------------

    %GT for orbit with no perturbations
    [gt_alpha,gt_delta,gt_long,gt_lat,gt_Y] =...
        gt_noPert(orb_y0,gt_longGreen,gt_tspan,Earth_mu,gt_we,gt_t0);

    %GT for orbit with J2 perturbation
    [gtJ2_alpha,gtJ2_delta,gtJ2_long,gtJ2_lat,gtJ2_Y] =...
        gt_J2(orb_y0,gt_longGreen,gt_tspan,Earth_mu,Earth_R,gt_we,gt_t0,j2);

    %Repeating GT for orbit with no perturbations
    m = 3; %rot of planet
    k = 4; % we have to choose them. %rev of sc
    
    gtRep_a              = gtRep_noPert(k,m,gt_we,Earth_mu);
    gtRep_kep0           = orb_kep0;
    gtRep_kep0(1)        = gtRep_a;
    gtRep_T              = 2*pi*sqrt(gtRep_a^3/Earth_mu);
    [gtRep_r0, gtRep_v0] = UT_kep2car(gtRep_kep0,Earth_mu);
    gtRep_y0             = [gtRep_r0; gtRep_v0]';
    gtRep_tspan          = linspace(0,gt_nOrb*gtRep_T,100000);
    
    [gtRep_alpha,gtRep_delta,gtRep_long,gtRep_lat] = ...
        gt_noPert(gtRep_y0,gt_longGreen,gtRep_tspan,Earth_mu,gt_we,gt_t0);

    %Repeating groundtrack with J2 perturbation
    gtRj2_a              =...
        gtRep_J2(k,m,gt_we,Earth_mu,j2,Earth_R,orb_kep0(1),orb_kep0(2),orb_kep0(3));
    gtRj2_kep0           = orb_kep0;
    gtRj2_kep0(1)        = gtRj2_a;
    [gtRj2_r0, gtRj2_v0] = UT_kep2car(gtRj2_kep0,Earth_mu);
    gtRj2_y0             = [gtRj2_r0; gtRj2_v0]';
    gtRj2_T              = 2*pi*sqrt(gtRj2_a^3/Earth_mu);
    gtRj2_tspan          = linspace(0,gt_nOrb*gtRj2_T,100000);
    
    [gtRj2_alpha,gtRj2_delta,gtRj2_long,gtRj2_lat] = ...
        gt_J2(gtRj2_y0,gt_longGreen,gtRj2_tspan,Earth_mu,Earth_R,gt_we,gt_t0,j2);

    % Conditions to avoid having horizontal lines connecting opposite longitudes
    for k = 2:length(gt_long)
        if gt_long(k)*gt_long(k-1)<0
            gt_long(k) = NaN;
        end
        if gtRep_long(k)*gtRep_long(k-1)<0
            gtRep_long(k) = NaN;
        end
        if gtJ2_long(k)*gtJ2_long(k-1)<0
            gtJ2_long(k) = NaN;
        end
        if gtRj2_long(k)*gtRj2_long(k-1)<0
            gtRj2_long(k) = NaN;
        end
    end

    %% Groundtrack plots

    %Groundtrack - no perturbation ------------------------------------------
    figure(1)
    c = imread('EarthTexture.jpg');
    c1 = imrotate(c,-180);
    image([180,-180],[-90,+90],c1)
    ax = gca;
    ax.YDir = 'normal';
    hold on
    plot(gt_long,gt_lat,'g')
    hold on 
    plot(gt_long(1),gt_lat(1),'go','Linewidth',2)
    hold on
    plot(gt_long(end),gt_lat(end),'gs','Linewidth',2)
    legend('Orbit','Start','End')
    title('Orbit without perturbation')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')

    %Groundtrack - j2 perturbation -------------------------------------------
    figure(2)
    c = imread('EarthTexture.jpg');
    c1 = imrotate(c,-180);
    image([180,-180],[-90,+90],c1)
    ax = gca;
    ax.YDir = 'normal';
    hold on
    plot(gt_long,gt_lat,'g')
    hold on 
    plot(gt_long(1),gt_lat(1),'go','Linewidth',2)
    hold on
    plot(gt_long(end),gt_lat(end),'gs','Linewidth',2)
    hold on
    plot(gtJ2_long,gtJ2_lat,'y')
    hold on
    plot(gtJ2_long(1),gtJ2_lat(1),'yo','Linewidth',2)
    hold on
    plot(gtJ2_long(end),gtJ2_lat(end),'ys','Linewidth',2)
    legend('Orbit','Start','End')
    title('Orbit with J2 perturbation')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')

    %Groundtrack - repeating - no perturbation -------------------------------
    figure(3)
    c = imread('EarthTexture.jpg');
    c1 = imrotate(c,-180);
    image([180,-180],[-90,+90],c1)
    ax = gca;
    ax.YDir = 'normal';
    hold on
    plot(gt_long,gt_lat,'g')
    hold on 
    plot(gt_long(1),gt_lat(1),'go','Linewidth',2)
    hold on
    plot(gt_long(end),gt_lat(end),'gs','Linewidth',2)
    hold on
    plot(gtRep_long,gtRep_lat,'r')
    hold on
    plot(gtRep_long(1),gtRep_lat(1),'ro','Linewidth',2)
    hold on
    plot(gtRep_long(end),gtRep_lat(end),'rs','Linewidth',2)
    legend('Orbit','Repeating orbit','Start','End','Start rep','End rep')
    title('Repeating groundtrack')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')

    %Groundtrack - repeating - j2 perturbation -------------------------------
    figure(4)
    c = imread('EarthTexture.jpg');
    c1 = imrotate(c,-180);
    image([180,-180],[-90,+90],c1)
    ax = gca;
    ax.YDir = 'normal';
    hold on
    plot(gt_long,gt_lat,'g')
    hold on 
    plot(gt_long(1),gt_lat(1),'go','Linewidth',2)
    hold on
    plot(gt_long(end),gt_lat(end),'gs','Linewidth',2)
    hold on
    plot(gtRep_long,gtRep_lat,'r')
    hold on 
    plot(gtRep_long(1),gtRep_lat(1),'ro','Linewidth',2)
    hold on
    plot(gtRep_long(end),gtRep_lat(end),'rs','Linewidth',2)
    hold on
    plot(gtRj2_long,gtRj2_lat,'y')
    hold on
    plot(gtRj2_long(1),gtRj2_lat(1),'yo','Linewidth',2)
    hold on
    plot(gtRj2_long(end),gtRj2_lat(end),'ys','Linewidth',2)
    legend('Original orbit','Start','End','Repeating Unperturbed','Start rep','End rep','Repeating perturbed','Start rep pert','End rep pert')
    title('Repeating groundtrack')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')


case 2 % Orbit propagations-----------------------------------------------
    %---Options for the ODE solver & Timespan
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
    
    orbProp_n = 3;
    switch orbProp_n
        case 1
            orb_n = 1;	    %number of orbits to simulate (1 period)
            orbProp_fil = 1;
        case 2
            orb_n = 36;	    %number of orbits to simulate (Moon period)
            orbProp_fil = 2;
        case 3
            orb_n = 500;	%number of orbits to simulate (1 years)
            orbProp_fil = 3;
        case 4
            orb_n = 18*500;	%number of orbits to simulate (18 years)
            orbProp_fil = 4;
    end
    
    tmax          = orb_n*orb_T;              %max time of simulation (nr of orbits*orbit period), in seconds.

    orb_stepsn    = 100;                      %number of segments per full orbit
    orb_stepst    = orb_T/orb_stepsn;         %time per step
    orb_stepsnTot = orb_stepsn*orb_n;         %total number of steps for the simulation
    orb_stepsYrs  = floor(tmax/(3600*24*365));

    tspan = linspace(0, tmax, orb_stepsnTot);

    %---Orbit propagation in cartesian coordinates (with J2+Moon perturbations)
    tic
    [orb_car_t, orb_car] = ode113(@(t,y)... 
        orbProp_carJ2Moon(t,y,Earth_mu,Moon_mu,j2,Earth_R,date0),...
        tspan, orb_y0, options);
    orb_car_r = [orb_car(:,1) orb_car(:,2) orb_car(:,3)];
    orb_car_v = [orb_car(:,4) orb_car(:,5) orb_car(:,6)];
    simTime_Cart = toc;

    % transformation to keplerian elements
     orb_car2kep = zeros(length(orb_car_t),6);
     for i=1:length(orb_car_t)
         orb_car2kep(i,:) = UT_car2kep(orb_car_r(i,:),orb_car_v(i,:),Earth_mu);
     end
     orb_car2kep(:,6) = unwrap(orb_car2kep(:,6));
     orb_car2kep(:,3:6) = rad2deg(orb_car2kep(:,3:6));

    %---Orbit propagation in keplerian elements
    tic
    [orb_kep_t, orb_kep] = ode113(@(t,y)...
        orbProp_gaussJ2Moon(t,y,Earth_mu,Moon_mu,j2,Earth_R,date0),...
        tspan,orb_kep0,options);
    simTime_Kep = toc;

    orb_kep(:,3:6) = rad2deg(orb_kep(:,3:6));

    %% Filtering--------------------------------------------------------------
    if orbProp_fil == 1 %n=1
        orb_kepFil = zeros(length(orb_kep_t),6);
        orb_kepFil(:,:) = NaN;
        
        orb_kepFil_long = zeros(length(orb_kep_t),6);
        orb_kepFil_long(:,:) = NaN;
        
    elseif orbProp_fil == 2 %n=36
        npoint_O   = nearest(3*orb_T/(sum(diff(orb_kep_t))/(numel(orb_kep_t)-1)));
        orb_kepFil = movmean(orb_kep,npoint_O,1); 
        %orb_kepFiltered(:,3:6) = rad2deg(orb_kepFiltered(:,3:6));
        
        orb_kepFil_long = zeros(length(orb_kep_t),6);
        orb_kepFil_long(:,:) = NaN;
        
    elseif orbProp_fil == 3 %n=500
        npoint_O   = nearest(3*orb_T/(sum(diff(orb_kep_t))/(numel(orb_kep_t)-1)));
        orb_kepFil = movmean(orb_kep,npoint_O,1); 
        %orb_kepFiltered(:,3:6) = rad2deg(orb_kepFiltered(:,3:6));
        
        orb_M           = 27*24*3600;
        npoint_M        = nearest(3*orb_M/(sum(diff(orb_kep_t))/(numel(orb_kep_t)-1)));
        orb_kepFil_long = movmean(orb_kep,npoint_M,1);
        
    elseif orbProp_fil == 4
        orb_kepFil = zeros(length(orb_kep_t),6);
        orb_kepFil(:,:) = NaN;%zeroes(6,length());
        
        orb_kepFil_long = zeros(length(orb_kep_t),6);
        orb_kepFil_long(:,:) = NaN;
    else
        orb_kepFil = NaN;
        
        orb_kepFil_long = NaN;
    end
    
    %% Method comparison------------------------------------------------------

    % Accuracy comparison
    error(:,1) = abs(orb_car2kep(:,1)-orb_kep(:,1))/orb_kep0(1);
    error(:,2) = abs(orb_car2kep(:,2)-orb_kep(:,2));
    error(:,3) = abs(orb_car2kep(:,3)-orb_kep(:,3))/(2*pi);
    error(:,4) = abs(orb_car2kep(:,4)-orb_kep(:,4))/(2*pi);
    error(:,5) = abs(orb_car2kep(:,5)-orb_kep(:,5))/(2*pi);
    error(:,6) = abs(orb_car2kep(:,6)-orb_kep(:,6))./abs(orb_kep(:,6));

    simAcc = max(error);

    %% Keplerian element evolution ploting
    figure(5)
    subplot(2,3,1);
    plot(tspan/orb_T,orb_kep(:,1),'b');
    hold on
    plot(tspan/orb_T,orb_car2kep(:,1),'r');
    hold on
    plot(tspan/orb_T,orb_kepFil(:,1),'g');
    hold on
    plot(tspan/orb_T,orb_kepFil_long(:,1),'m');
    title('semi-major axis')
    xlabel('Time [T]')
    ylabel('a [km]')
    if orbProp_n == 1 || orbProp_n == 4
        legend ('Gauss equations','Cartesian method')
    elseif orbProp_n == 2
        legend ('Gauss equations','Cartesian method','Long period')
    elseif orbProp_n == 3
        legend ('Gauss equations','Cartesian method','Long period','Secular')
    end

    subplot(2,3,2);
    plot(tspan/orb_T,orb_kep(:,2),'b');
    hold on
    plot(tspan/orb_T,orb_car2kep(:,2),'r');
    hold on
    plot(tspan/orb_T,orb_kepFil(:,2),'g');
    hold on
    plot(tspan/orb_T,orb_kepFil_long(:,2),'m');
    title('Excentricity');
    xlabel('Time [T]')
    ylabel('e [-]')
    if orbProp_n == 1 || orbProp_n == 4
        legend ('Gauss equations','Cartesian method')
    elseif orbProp_n == 2
        legend ('Gauss equations','Cartesian method','Secular')
    elseif orbProp_n == 3
        legend ('Gauss equations','Cartesian method','Long period','Secular')
    end

    subplot(2,3,3);
    plot(tspan/orb_T,orb_kep(:,3),'b');
    hold on
    plot(tspan/orb_T,orb_car2kep(:,3),'r');
    hold on
    plot(tspan/orb_T,orb_kepFil(:,3),'g');
    hold on
    plot(tspan/orb_T,orb_kepFil_long(:,3),'m');
    title('inclination')
    xlabel('Time [T]')
    ylabel('i [º]')
    if orbProp_n == 1 || orbProp_n == 4
        legend ('Gauss equations','Cartesian method')
    elseif orbProp_n == 2
        legend ('Gauss equations','Cartesian method','Secular')
    elseif orbProp_n == 3
        legend ('Gauss equations','Cartesian method','Long period','Secular')
    end

    subplot(2,3,4);
    plot(tspan/orb_T,orb_kep(:,4),'b');
    hold on
    plot(tspan/orb_T,orb_car2kep(:,4),'r');
    hold on
    plot(tspan/orb_T,orb_kepFil(:,4),'g');
    hold on
    plot(tspan/orb_T,orb_kepFil_long(:,4),'m');
    title('Right Ascension of Ascending Node')
    xlabel('Time [T]')
    ylabel('RAAN [º]')
    if orbProp_n == 1 || orbProp_n == 4
        legend ('Gauss equations','Cartesian method')
    elseif orbProp_n == 2
        legend ('Gauss equations','Cartesian method','Secular')
    elseif orbProp_n == 3
        legend ('Gauss equations','Cartesian method','Long period','Secular')
    end

    subplot(2,3,5);
    plot(tspan/orb_T,orb_kep(:,5),'b');
    hold on
    plot(tspan/orb_T,orb_car2kep(:,5),'r');
    hold on
    plot(tspan/orb_T,orb_kepFil(:,5),'g');
    hold on
    plot(tspan/orb_T,orb_kepFil_long(:,5),'m');
    title('Argument of perigee')
    xlabel('Time [T]')
    ylabel('w [º]')
    if orbProp_n == 1 || orbProp_n == 4
        legend ('Gauss equations','Cartesian method')
    elseif orbProp_n == 2
        legend ('Gauss equations','Cartesian method','Secular')
    elseif orbProp_n == 3
        legend ('Gauss equations','Cartesian method','Long period','Secular')
    end

    subplot(2,3,6);
    plot(tspan/orb_T,orb_kep(:,6),'b');
    hold on
    plot(tspan/orb_T,orb_car2kep(:,6),'r');
    hold on
    plot(tspan/orb_T,orb_kepFil(:,6),'g');
    hold on
    plot(tspan/orb_T,orb_kepFil_long(:,6),'m');
    title('true anomaly')
    xlabel('Time [T]')
    ylabel('f [º]')
    if orbProp_n == 1 || orbProp_n == 4
        legend ('Gauss equations','Cartesian method')
    elseif orbProp_n == 2
        legend ('Gauss equations','Cartesian method','Secular')
    elseif orbProp_n == 3
        legend ('Gauss equations','Cartesian method','Long period','Secular')
    end

    %% Error plotting

    figure(6)
    subplot(2,3,1);
    semilogy(tspan/orb_T,error(:,1),'b');
    grid on
    title('semi-major axis')
    xlabel('Time [T]')
    ylabel('a [km]')

    subplot(2,3,2);
    semilogy(tspan/orb_T,error(:,2),'b');
    grid on
    title('Excentricity');
    xlabel('Time [T]')
    ylabel('e [-]')

    subplot(2,3,3);
    semilogy(tspan/orb_T,error(:,3),'b');
    grid on
    title('inclination')
    xlabel('Time [T]')
    ylabel('i [º]')

    subplot(2,3,4);
    semilogy(tspan/orb_T,error(:,4),'b');
    grid on
    title('RAAN')
    xlabel('Time [T]')
    ylabel('RAAN [º]')

    subplot(2,3,5);
    semilogy(tspan/orb_T,error(:,5),'b');
    grid on
    title('Argument of perigee')
    xlabel('Time [T]')
    ylabel('w [º]')

    subplot(2,3,6);
    semilogy(tspan/orb_T,error(:,6),'b');
    grid on
    title('true anomaly')
    xlabel('Time [T]')
    ylabel('f [º]')


case 3 % Orbit evolution plot
    UT_orbPlot


case 4 % Comparison with real data    
    
    %Reading from the file & Separating keplerian elements----------------
    rd_source = 1; %'1' = data for 1 year 
                   %'2' = data for 1 day
    switch rd_source
        case 1
            rd_data = csvread('horizons-results-long.txt',0,2);
            rd_kep = rd_data(:,[10 1 3 4 5 9]);
            rd_time = 3600*24;
            rd_tmax = 3600*24*365;

        case 2
            rd_data = csvread('horizons-results-short.txt',0,2);
            rd_kep = rd_data(:,[10 1 3 4 5 9]);
            rd_time = 3600;
            rd_tmax = 3600*24;

    end
    
    rd_tspan = linspace(0,rd_tmax,length(rd_kep));
    
    rd_kep0 = rd_kep(1,:);                                                 
    rd_kep0(3:6) = deg2rad(rd_kep0(3:6));                                  
    %rd_kep(:,6) = unwrap(rd_kep0(:,6));                                         

    %Propagation from initial conditions ---------------------------------
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);                    
    rdProp_stepsn = 100;                                                   
    rdProp_tspan  = linspace(0,rd_tmax,10000);%rdProp_stepsn*(length(rd_kep)-1));  

    [rd_kep_t, rdProp_kep] = ode113(@(t,y)...                              
        orbProp_gaussJ2Moon(t,y,Earth_mu,Moon_mu,j2,Earth_R,date0),...
        rdProp_tspan,rd_kep0,options);

    rdProp_kep(:,3:6) = rad2deg(rdProp_kep(:,3:6));                        
    rdProp_kep(:,6) = wrapTo360(rdProp_kep(:,6));                          

    %Ploting--------------------------------------------------------------
    figure(8)
    subplot(2,3,1);
    plot(rd_tspan/rd_time,rd_kep(:,1),'o');
    hold on
    plot(rdProp_tspan/rd_time,rdProp_kep(:,1),'r');
    title('semi-major axis')
    if rd_source == 1
        xlabel('days')
    elseif rd_source == 2
        xlabel('hours')
    end
    ylabel('a [km]')
    legend ('Real data','Gauss equations')

    subplot(2,3,2);
    plot(rd_tspan/rd_time,rd_kep(:,2),'o');
    hold on
    plot(rdProp_tspan/rd_time,rdProp_kep(:,2),'r');
    title('excentricity')
    if rd_source == 1
        xlabel('days')
    elseif rd_source == 2
        xlabel('hours')
    end
    ylabel('[-]')

    subplot(2,3,3);
    plot(rd_tspan/rd_time,rd_kep(:,3),'o');
    hold on
    plot(rdProp_tspan/rd_time,rdProp_kep(:,3),'r');
    title('inclination')
    if rd_source == 1
        xlabel('days')
    elseif rd_source == 2
        xlabel('hours')
    end
    ylabel('i [º]')

    subplot(2,3,4);
    plot(rd_tspan/rd_time,rd_kep(:,4),'o');
    hold on
    plot(rdProp_tspan/rd_time,rdProp_kep(:,4),'r');
    title('RAAN')
    if rd_source == 1
        xlabel('days')
    elseif rd_source == 2
        xlabel('hours')
    end
    ylabel('RAAN [º]')

    subplot(2,3,5);
    plot(rd_tspan/rd_time,rd_kep(:,5),'o');
    hold on
    plot(rdProp_tspan/rd_time,rdProp_kep(:,5),'r');
    title('w')
    if rd_source == 1
        xlabel('days')
    elseif rd_source == 2
        xlabel('hours')
    end
    ylabel('w [º]')

    subplot(2,3,6);
    plot(rd_tspan/rd_time,rd_kep(:,6),'o');
    hold on
    plot(rdProp_tspan/rd_time,rdProp_kep(:,6),'r');
    title('True anomaly')
    if rd_source == 1
        xlabel('days')
    elseif rd_source == 2
        xlabel('hours')
    end
    ylabel('f [º]')

end