function [r,v] =UT_kep2car(kep_coord,mu_e)

%kep2car computes position and velocity in keplerian frame as vectors starting from keplerian parameters.

%PROTOTYPE: 
%     [r,v] =kep2car(kepl_coord,mu)
% 
% INPUT:
%     kep_coord [6x1]     keplerian coordinates = [a,e,i,omega,w,f], [km,~,rad,rad,rad,rad]
%     mu [1]              Gravitational constant of the Earth [km^3/s^2]
%     
% OUTPUT:
%     r [3x1]             Position vector [km]
%     v [3x1]             Velocity vectors [km/s]
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


a     = kep_coord(1);		%semi-major axis
e     = kep_coord(2);		%eccentricity
i     = kep_coord(3);		%inclination
omega = kep_coord(4);       %right ascension of ascending node
w     = kep_coord(5);		%argument of perigee
f     = kep_coord(6);		%true anomaly

h      = sqrt(a*(1-e^2)*mu_e);                              %angular momentum
r_peri = h^2/(mu_e*(1+e*cos(f)))*[cos(f); sin(f); 0];		%position vector in perifocal RF
v_peri = mu_e/h*[-sin(f); e+cos(f); 0];                     %velocity vector in perifocal RF

R3_w     = [cos(w), sin(w), 0; -sin(w), cos(w),0; 0,0,1];
R1_i     = [1,0,0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R3_omega = [cos(omega), sin(omega),0; -sin(omega), cos(omega), 0; 0,0,1];

r = (R3_w*R1_i*R3_omega)'*r_peri;
v = (R3_w*R1_i*R3_omega)'*v_peri;


