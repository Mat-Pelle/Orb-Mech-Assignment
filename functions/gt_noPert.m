function [gt_alpha,gt_delta,gt_lon,gt_lat,gt_Y] = gt_noPert(y0,gt_longGreen,gt_tspan,Earth_mu,gt_we,gt_t0)
%Groundtracks computes RAAN, declination, latitude and longitude to have the projection of 
%the orbit of s/c onto Earth's surface

%PROTOTYPE: 
%     [alpha,delta,lon,lat] = Groundtracks(kep_coord,long_green,t,mu,we,t0)
% 
% INPUT:
%     orb_y0 [6x1]        keplerian elements ([a,e,i,omega,w,f]), 
%                         [km,~,rad,rad,rad,rad]
%     gt_longGreen [1]    longitude of Greenwich meridian [rad]
%     gt_tspan            Time span of the orbit propagation[s]
%     Earth_mu [1]        Gravitational constant of the Earth [km^3/s^2]
%     gt_we [1]           Angular velocity of the Earth [rad/s]
%     gt_t0 [1]           Initial time [s]
%     
% OUTPUT:
%     gt_alpha               Right Ascension [deg]     
%     gt_delta               Declination [deg]
%     gt_lon                 Longitude [deg]
%     gt_lat                 Latitude [deg]
%     gt_Y                   Position vector of the spacecraft, in
%                            cartesian coordinates [km km km]
%
% CONTRIBUTORS
%       Bertolini Edoardo
%       Busi Silvia
%       Muylle Julia
%       Pellegrini Matias
%
% VERSIONS
%
% 30/11/2020: First Version
% 20/1/2021: Definitive version

options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[~,gt_Y]   = ode113(@(t,y) orbProp_car(t,y,Earth_mu),gt_tspan,y0,options);
gt_Y       = [gt_Y(:,1) gt_Y(:,2) gt_Y(:,3)];

%Right ascension of the ascending node and declination
for i =1:length(gt_Y(:,1))
    
    gt_delta(i) = asin(gt_Y(i,3)/norm(gt_Y(i,:)));
    gt_alpha(i) = acos(gt_Y(i,1)/(norm(gt_Y(i,:))*cos(gt_delta(i))));
    
    if (gt_Y(i,2)/norm(gt_Y(i,:))) <= 0
        gt_alpha(i) = 2*pi-gt_alpha(i);
    end
    
    %Longitude of Greenwich meridian
    thetag(i) = gt_longGreen+gt_we*(gt_tspan(i)-gt_t0);
    
    %Longitude
    gt_lon(i)    = gt_alpha(i)-thetag(i);
    gt_lon(i)    = wrapToPi(gt_lon(i));
    
    %Latitude
    gt_lat(i)    = gt_delta(i);
    
end

gt_lat = rad2deg(gt_lat);
gt_lon = rad2deg(gt_lon);

