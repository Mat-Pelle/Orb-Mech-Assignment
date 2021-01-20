function [gtJ2_alpha,gtJ2_delta,gtJ2_lon,gtJ2_lat,gtJ2_Y] = gt_J2(y0,long_green,tspan,Earth_mu,Earth_R,we,t0,j2)

% Groundtracks_J2 computes RAAN, declination, latitude and longitude to have the projection 
%of the orbit of s/c onto Earth’s surface considering J2 effect

%PROTOTYPE: 
%     [alpha,delta,lon,lat] = Groundtracks_J2(kep_coord,long_green,t,mu,we,t0)
% 
% INPUT:
%     orb_y0 [6x1]        keplerian elements ([a,e,i,omega,w,f]), 
%                         [km,~,rad,rad,rad,rad]
%     gt_longGreen [1]    longitude of Greenwich meridian [rad]
%     gt_tspan            Time span of the orbit propagation[s]
%     Earth_mu [1]        Gravitational constant of the Earth [km^3/s^2]
%     Earth_R  [1]        Radius of the Earth [km]
%     gt_we [1]           Angular velocity of the Earth [rad/s]
%     gt_t0 [1]           Initial time [s]
%     J2 [1]              J2 effect
%     
% OUTPUT:
%     gtJ2_alpha               Right Ascension [deg]     
%     gtJ2_delta               Declination [deg]
%     gtJ2_lon                 Longitude [deg]
%     gtJ2_lat                 Latitude [deg]
%     gtJ2_Y                   Position vector of the spacecraft, in
%                              cartesian coordinates [km km km]
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


%[r0,v0] = kep2car(kep_coord,mu);
%y0      = [r0,v0];
options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[~,gtJ2_Y]   = ode113(@(t,y) orbProp_carJ2(t,y,Earth_mu,Earth_R,j2),tspan,y0,options);
gtJ2_Y       = [gtJ2_Y(:,1) gtJ2_Y(:,2) gtJ2_Y(:,3)];

%Right ascension of the ascending node and declination
for i =1:length(gtJ2_Y(:,1))
    
    gtJ2_delta(i) = asin(gtJ2_Y(i,3)/norm(gtJ2_Y(i,:)));
    gtJ2_alpha(i) = acos(gtJ2_Y(i,1)/(norm(gtJ2_Y(i,:))*cos(gtJ2_delta(i))));
    
    if (gtJ2_Y(i,2)/norm(gtJ2_Y(i,:)))<=0
        gtJ2_alpha(i) = 2*pi-gtJ2_alpha(i);
    end
    %Longitude of Greenwich meridian
    thetag(i) = long_green+we*(tspan(i)-t0);
    
    %Longitude
    gtJ2_lon(i)    = gtJ2_alpha(i)-thetag(i);
    
    %Latitude
    gtJ2_lat(i)    = gtJ2_delta(i);
    gtJ2_lon(i)    = wrapToPi(gtJ2_lon(i));
end

gtJ2_lat = rad2deg(gtJ2_lat);
gtJ2_lon = rad2deg(gtJ2_lon);
