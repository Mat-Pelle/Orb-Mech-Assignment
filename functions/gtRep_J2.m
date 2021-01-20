function gtRj2_a = gtRep_J2(k,m,gt_we,Earth_mu,J2,Earth_R,a0,e,i)
% % rep_groundtrack_J2 computes the value of the semi-major axis to obtain repeating 
%groundtrack with J2 effect

%PROTOTYPE: 
%     a = rep_groundtrack_J2(k,m,we,mu,J2,Re,a0,e,i)
%
% INPUT:
%     k [1]          Revolution of the satelite
%     m [1]          Rotations of the planet
%     gt_we [1]      Angular velocity of the Earth [rad/s]
%     Earth_mu [1]   Gravitational constant of the Earth [km^3/s^2]
%     J2 [1]         J2 effect
%     a0 [1]         Initial semimajor axis [km]
%     e [1]          eccentricity
%     i [1]          inclination [rad]
%     
% OUTPUT:
%    a [1]        Semimajor axis [km]
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

%Secular evolution
n      = @(a) sqrt(Earth_mu/a^3); 
domega = @(a) -(3/2*sqrt(Earth_mu)*J2*Earth_R^2/((1-e^2)^2*a^(7/2)))*cos(i);
dw     = @(a) -(3/2*sqrt(Earth_mu)*J2*Earth_R^2/((1-e^2)^2*a^(7/2)))*(5/2*(sin(i))^2-2);
dM0    = @(a) (3/2*sqrt(Earth_mu)*J2*Earth_R^2/((1-e^2)^1.5*a^(7/2)))*(1-3/2*(sin(i))^2);

fun     = @(a) (gt_we-domega(a))/(n(a)+dw(a)+dM0(a))-m/k;
options = optimset('TolFun',1e-16,'Display','off');
gtRj2_a       = fsolve(fun, a0, options);


