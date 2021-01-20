function orbDot = orbProp_gaussJ2Moon(t, y, mu1, mu2, j2, Re, date0)

%INPUT:
% 
% OUTPUT:
% 

%% Start-up

a    = y(1);
e    = y(2);
i    = y(3);
RAAN = y(4);
w    = y(5);
f    = y(6);

[car_r, car_v] = UT_kep2car([a,e,i,RAAN,w,f],mu1);
rmod = norm(car_r);
vmod = norm(car_v);

%% Change of reference frame

%--- TNH
tt = car_v/vmod;
hh = cross(car_r,car_v);%/norm(cross(r,v));
hh = hh/norm(hh);
nn = cross(hh,tt);

%--- RSW
rr = car_r/norm(car_r);
ww = cross(car_r,car_v);%/norm(cross(r,v));
ww = ww/norm(ww);
ss = cross(ww,rr);
%ss = ss/norm(ss);

% car2tnh = [tt(:)'; nn(:)'; hh(:)']';

% tnh2rsw = [e*sin(f)    -(1+e*cos(f)) 0
%            1+e*cos(f)  e*sin(f)      0         
%            0           0             1];

% car2rsw = [rr(:)'; ss(:)'; ww(:)']';

%% Perturbations

b = a*sqrt(1-(e^2));
p = (b^2)/a;
r = p/(1+(e*cos(f)));
v = sqrt(((2*mu1)/r)-(mu1/a));
n = sqrt(mu1/(a^3));
h = n*a*b;

%--- J2 Perturbation
aj2 = (3/2)*((j2*mu1*Re^2)/rmod^4)*[((car_r(1)/rmod)*(5*(y(3)^2/rmod^2)-1))
                                    ((car_r(2)/rmod)*(5*(y(3)^2/rmod^2)-1))
                                    ((car_r(3)/rmod)*(5*(y(3)^2/rmod^2)-3))];

% aj2_tnh = car2tnh*aj2;                                
% aj2_rsw2 = tnh2rsw*aj2_tnh

% aj2_rsw3 = car2rsw*aj2
aj2_rsw4 = [dot(aj2,rr); dot(aj2,ss); dot(aj2,ww)];

aj2_rsw = -(3/2)*((j2*mu1*Re^2)/r^4)*[1-(3*(sin(i)^2)*((sin(f+w))^2))
                                         sin(i)^2*sin(2*(f+w))
                                         sin(2*i)*sin(f+w)];  

rsw_rr = aj2_rsw4(1)-aj2_rsw(1);
rsw_ss = aj2_rsw4(2)-aj2_rsw(2);                                     
rsw_ww = aj2_rsw4(3)-aj2_rsw(3);
rsw_error = [rsw_rr; rsw_ss; rsw_ww];

%--- Moon perturbation
date = date0+(t/(3600*24));
[Moon_r, ~] = ephMoon(date); %Initial Moon position and velocity vectors, in [km] & [km/s], respectively

%aMoon = mu2*(((y-Moon_r)/((norm(r-Moon_r))^3))-(Moon_r/((norm(Moon_r))^3)));
 aMoon2 = mu2*[(((car_r(1)-Moon_r(1))/((norm(car_r-Moon_r))^3))-(Moon_r(1)/((norm(Moon_r))^3)))
              (((car_r(2)-Moon_r(2))/((norm(car_r-Moon_r))^3))-(Moon_r(2)/((norm(Moon_r))^3)))
              (((car_r(3)-Moon_r(3))/((norm(car_r-Moon_r))^3))-(Moon_r(3)/((norm(Moon_r))^3)))];
          
%alternative calculation of the moon perturbation
Moon_rel = Moon_r-car_r';
%Moon_rel = Moon_rel/norm(Moon_rel);
q = dot(car_r,((2*Moon_r - car_r')/(norm(Moon_r))^2));
F = (q^2 - 3*q + 3)*q/(1 + (1-q)^1.5);

aMoon = mu2/(norm(Moon_rel))^3*(F*Moon_r - car_r');

%ax = aj2(1)+aMoon(1);
%ay = aj2(2)+aMoon(2);
%az = aj2(3)+aMoon(3);
%a_car = [ax ay az];
a_car = aj2+aMoon'

ar = dot(a_car,rr);
as = dot(a_car,ss);
aw = dot(a_car,ww);

% a_tnh = car2tnh*a_xyz;
% a_rsw = tnh2rsw*a_tnh;

% ar = a_rsw(1);
% as = a_rsw(2);
% aw = a_rsw(3);

% aMoon_tnh = car2tnh'*aMoon;

% aMoon_rsw = tnh2rsw'*aMoon_tnh;          

%--- Accelerations due to perturbations
% at = aj2_tnh(1)+aMoon_tnh(1);
% an = aj2_tnh(2)+aMoon_tnh(2);
% ah = aj2_tnh(3)+aMoon_tnh(3);

% ar = aj2_rsw(1)+aMoon_rsw(1);
% as = aj2_rsw(2)+aMoon_rsw(2);
% aw = aj2_rsw(3)+aMoon_rsw(3);

%% Keplerian element variations

%--- in TNH
% aDot     = ((2*(a^2)*v)/Earth_mu)*at;
% eDot     = (1/v)*((2*(e+cos(f))*at)-((r/a)*sin(f*an)));
% iDot     = ((r*cos(f+w))/h)*ah;
% OmegaDot = ((r*sin(f+w))/(h*sin(i)))*ah;
% wDot     = ((1/(e*v))*((2*sin(f*at))+(an*((2*e)+((r/a)*cos(f))))))-(ah*((r*sin(f+w)*cos(i))/(h*sin(i))));
% fDot     = (h/(r^2))-((1/(e*v))*((2*sin(f*at))+(an*((2*e)+((r/a)*cos(f))))));
% % MDot     = n-((b/(e*a*v))*((2*(1+((exp(2*r))/p))*sin(f*at))+((r/a)*cos(f*an)))); 

%--- in RSW
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