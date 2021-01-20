function orbDot = orbProp_gaussJ2Moon(t, y, Earth_mu, Moon_mu, j2, Earth_R, date0)

%orbProp_carJ2Moon provides the derivative of the state vector using the 
%gaussian equations and taking into account J2 and Moon perturbations.
%
%PROTOTYPE: 
%   dy = odefun_J2 (t,y,mu,J2)
% 
%INPUT:
%   t        [1]   Time [s]
%   y        [6x1] State Vector containing the keplerian elements
%   Earth_mu [1]   Gravitational constant of the Earth [km^3/s^2]
%   Moon_mu  [1]   Gravitational constant of the Moon [km^3/s^2]
%   J2       [1]   J2 effect
%   Earth_R  [1]   Radius of the Earth [km]
%   date0    [1]   Starting date of the simulation
%     
% OUTPUT:
%   orbDot [6x1]     Derivative fo the state vector
%
% CONTRIBUTORS
%   Bertolini Edoardo
%   Busi Silvia
%   Muylle Julia
%   Pellegrini Matias
%
% VERSIONS
%   18/01/2021: Definitive Version

%% Start-up

a    = y(1);
e    = y(2);
i    = y(3);
RAAN = y(4);
w    = y(5);
f    = y(6);

[car_r, car_v] = UT_kep2car([a,e,i,RAAN,w,f],Earth_mu);
car_r = car_r';
rmod = norm(car_r);
vmod = norm(car_v);

%% RSW reference frame

rr = car_r;%/norm(car_r);
rr = rr/norm(rr);
ww = cross(car_r,car_v);%/norm(cross(r,v));
ww = ww/norm(ww);
ss = cross(ww,rr);
%ss = ss/norm(ss);

%% Perturbations

b = a*sqrt(1-(e^2));
p = (b^2)/a;
r = p/(1+(e*cos(f)));
v = sqrt(((2*Earth_mu)/r)-(Earth_mu/a));
n = sqrt(Earth_mu/(a^3));
h = n*a*b;

%--- J2 Perturbation
aj2 = (3/2)*((j2*Earth_mu*Earth_R^2)/r^4)*[((car_r(1)/r)*(5*(car_r(3)^2/r^2)-1))
                                      ((car_r(2)/r)*(5*(car_r(3)^2/r^2)-1))
                                      ((car_r(3)/r)*(5*(car_r(3)^2/r^2)-3))];
                                                         
aj2_rsw = [dot(aj2,rr); dot(aj2,ss); dot(aj2,ww)];

% aj2_rsw = -(3/2)*((j2*muE*Re^2)/r^4)*[1-(3*((sin(i))^2)*((sin(f+w))^2))
%                                          sin(i)^2*sin(2*(f+w))
%                                          sin(2*i)*sin(f+w)];  

% rsw_rr = aj2_rsw4(1)-aj2_rsw(1);
% rsw_ss = aj2_rsw4(2)-aj2_rsw(2);                                     
% rsw_ww = aj2_rsw4(3)-aj2_rsw(3);
% rsw_error = [rsw_rr; rsw_ss; rsw_ww]

%--- Moon perturbation
date = date0+(t/(3600*24));
[Moon_r, ~] = ephMoon(date); %Initial Moon position and velocity vectors, in [km] & [km/s], respectively

aMoon = Moon_mu*(((car_r-Moon_r)/((norm(car_r-Moon_r))^3))-(Moon_r/((norm(Moon_r))^3)));

% aMoon2 = muM*[(((car_r(1)-Moon_r(1))/((norm(car_r-Moon_r))^3))-(Moon_r(1)/((norm(Moon_r))^3)))
%               (((car_r(2)-Moon_r(2))/((norm(car_r-Moon_r))^3))-(Moon_r(2)/((norm(Moon_r))^3)))
%               (((car_r(3)-Moon_r(3))/((norm(car_r-Moon_r))^3))-(Moon_r(3)/((norm(Moon_r))^3)))];
% 
% aMoon_err = aMoon-aMoon2'
          
%alternative calculation of the moon perturbation------------------------
% Moon_rel = Moon_r-car_r;
% Moon_rel = Moon_rel/norm(Moon_rel);
% q = dot(car_r,((2*Moon_r - car_r)/(norm(Moon_r))^2));
% F = (q^2 - 3*q + 3)*q/(1 + (1-q)^1.5);

%aMoon = mu2/(norm(Moon_rel))^3*(F*Moon_r - car_r');

aMoon_rsw(1) = dot(aMoon,rr);
aMoon_rsw(2) = dot(aMoon,ss);
aMoon_rsw(3) = dot(aMoon,ww);

ar = aj2_rsw(1)+aMoon_rsw(1);
as = aj2_rsw(2)+aMoon_rsw(2);
aw = aj2_rsw(3)+aMoon_rsw(3);

%% Keplerian element variations

aDot     = ((2*a^2)/h)*((e*sin(f)*ar)+(p*as/r));
% hDot     = r*as;
eDot     = (1/h)*((p*sin(f)*ar)+((((p+r)*cos(f))+r*e)*as));
iDot     = (r*cos(f+w))*aw/h;
OmegaDot = (r*sin(f+w))*aw/(h*sin(i));
wDot     = (1/(h*e))*((-p*cos(f)*ar)+((p+r)*sin(f)*as))-(((r*sin(f+w)*cos(i))*aw/(h*sin(i))));
fDot     = (h/(r^2))+((1/(h*e))*((p*cos(f)*ar)-((p+r)*sin(f)*as)));
% MDot     = n + ((b/(a*h*e))*(ar*((p*cos(f))-(2*r*e))-((p+r)*(sin(f*as)))));

orbDot = [aDot
          eDot
          iDot
          OmegaDot
          wDot
          fDot];

end