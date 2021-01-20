function dy = orbProp_carJ2Moon(t,y,Earth_mu,Moon_mu,j2,Earth_R,date0)

%orbProp_carJ2Moon provides the derivative of the state vector using the 
%dynamic equation of 2BP with J2 and Moon perturbation
%
%PROTOTYPE: 
%   dy = odefun_J2 (t,y,mu,J2)
% 
%INPUT:
%   t              Time [s]
%   y        [6x1] State Vector containing position and velocity vectors concatanate
%   Earth_mu [1]   Gravitational constant of the Earth [km^3/s^2]
%   Moon_mu  [1]   Gravitational constant of the Moon [km^3/s^2]
%   J2       [1]   J2 effect
%   Earth_R  [1]   Radius of the Earth [km]
%   date0    [1]   Starting date of the simulation
%     
% OUTPUT:
%   dy [6x1]     Derivative fo the state vector
%
% CONTRIBUTORS
%   Bertolini Edoardo
%   Busi Silvia
%   Muylle Julia
%   Pellegrini Matias
%
% VERSIONS
%   18/01/2021: Definitive Version

r      = y(1:3);
r_norm = norm(r);
%Re     = UT_astroConstants(23);
k      = 3/2*j2*Earth_mu*Earth_R^2/r_norm^4;

aj(1) = k*(r(1)./r_norm.*(5*r(3).^2./r_norm^2-1));
aj(2) = k*(r(2)./r_norm.*(5*r(3).^2./r_norm^2-1));
aj(3) = k*(r(3)./r_norm.*(5*r(3).^2./r_norm^2-3));

date0 = date0+(t/(3600*24));
[Moon_r, ~] = ephMoon(date0); %Initial Moon position and velocity vectors, in [km] & [km/s], respectively

aMoon = Moon_mu*(((Moon_r'-r)/((norm(Moon_r'-r))^3))-(Moon_r'/((norm(Moon_r'))^3)));

% aMoon(1) = mu2*(((r(1)-Moon_r(1))/((norm(r-Moon_r))^3))-(Moon_r(1)/((norm(Moon_r))^3)));
% aMoon(2) = mu2*(((r(2)-Moon_r(2))/((norm(r-Moon_r))^3))-(Moon_r(2)/((norm(Moon_r))^3)));
% aMoon(3) = mu2*(((r(3)-Moon_r(3))/((norm(r-Moon_r))^3))-(Moon_r(3)/((norm(Moon_r))^3)));



dy = [y(4:6); -Earth_mu/(norm(r)^3).*y(1:3)+ aj'+aMoon];
end