function orbDot = orbProp_gauss1(y, mu)

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

[car_r car_v] = UT_kep2car([a,e,i,RAAN,w,f],mu);

%car_r = [car(1) car(2) car(3)];
%car_v = [car(4) car(5) car(6)];

rmod = norm(car_r);
vmod = norm(car_v);

%% Change of reference frame

% %TNH
% tt = car_v/vmod;
% hh = cross(car_r,car_v);%/norm(cross(r,v));
% hh = hh/norm(hh);
% nn = cross(hh,tt);
% 
% %RSW
% rr = car_r/rmod;
% ww = cross(car_r,car_v);%/norm(cross(r,v));
% ww = ww/norm(ww);
% ss = cross(ww,rr);
% 
% car2tnh = [tt; nn; hh];
% 
% tnh2rsw = [e*sin(f)    -(1+e*cos(f)) 0
%            1+e*cos(f)  e*sin(f)      0         
%            0           0             1];

%% Perturbations

b = a*sqrt(1-(e^2));
p = (b^2)/a;
r = p/(1+(e*cos(f)));
v = sqrt(((2*mu)/r)-(mu/a));
n = sqrt(mu/(a^3));
h = sqrt(p*mu);%n*a*b;

%--- J2 Perturbation
% aj2 = (3/2)*((j2*mu*re^2)/r_mod^4)*[((car_r(1)/r_mod)*(5*(y(3)^2/r_mod^2)-1))
%                                     ((car_r(2)/r_mod)*(5*(y(3)^2/r_mod^2)-1))
%                                     ((car_r(3)/r_mod)*(5*(y(3)^2/r_mod^2)-3))];

%aj2_tnh = car2tnh*aj2;                              
                                
% aj2_rsw = -(3/2)*((j2*Earth_mu*Earth_R^2)/r^4)*[1-(3*sin(i)^2*sin(f+w)^2)
%                                                    sin(i)^2*sin(2*(f+w))
%                                                    sin(2*i)*sin(f+w)];  

%Moon
% date = (date0/(3600*24))+t;
% [Moon_r, ~] = ephMoon(date); %Initial Moon position and velocity vectors, in [km] & [km/s], respectively
% 
% aMoon = Moon_mu*[(((car_r(1)-Moon_r(1))/((norm(car_r-Moon_r))^3))-(Moon_r(1)/((norm(Moon_r))^3)))
%                   (((car_r(2)-Moon_r(2))/((norm(car_r-Moon_r))^3))-(Moon_r(2)/((norm(Moon_r))^3)))
%                   (((car_r(3)-Moon_r(3))/((norm(car_r-Moon_r))^3))-(Moon_r(3)/((norm(Moon_r))^3)))];
%               
% aMoon_tnh = car2tnh*aMoon;
% 
% aMoon_rsw = tnh2rsw*aMoon_tnh;          

%Accelerations due to perturbations
% at = aj2_tnh(1)+aMoon_tnh(1);
% an = aj2_tnh(2)+aMoon_tnh(2);
% ah = aj2_tnh(3)+aMoon_tnh(3);

ar = 0;%aj2_rsw(1);%+aMoon_rsw(1);
as = 0;%aj2_rsw(2);%+aMoon_rsw(2);
aw = 0;%aj2_rsw(3);%+aMoon_rsw(3);

%% Main function

%in TNH
% aDot     = ((2*(a^2)*v)/Earth_mu)*at;
% eDot     = (1/v)*((2*(e+cos(f))*at)-((r/a)*sin(f*an)));
% iDot     = ((r*cos(f+w))/h)*ah;
% OmegaDot = ((r*sin(f+w))/(h*sin(i)))*ah;
% wDot     = ((1/(e*v))*((2*sin(f*at))+(an*((2*e)+((r/a)*cos(f))))))-(ah*((r*sin(f+w)*cos(i))/(h*sin(i))));
% fDot     = (h/(r^2))-((1/(e*v))*((2*sin(f*at))+(an*((2*e)+((r/a)*cos(f))))));
% % MDot     = n-((b/(e*a*v))*((2*(1+((exp(2*r))/p))*sin(f*at))+((r/a)*cos(f*an)))); 

%in RSW
aDot     = ((2*a^2)/h)*((e*sin(f)*ar)+((p/r)*as));
% hDot     = r*as;
eDot     = (1/h)*((p*sin(f)*ar)+((((p+r)*cos(f))+r*e)*as));
iDot     = ((r*cos(f+w))/h)*aw;
OmegaDot = ((r*sin(f+w))/(h*sin(i)))*aw;
wDot     = (1/(h*e))*((-p*cos(f)*ar)+((p+r)*sin(f)*as))-(((r*sin(f+w)*cos(i))/(h*sin(i)))*aw);
fDot     = (h/(r^2))+((1/(h*e))*((p*cos(f)*ar)-((p+r)*sin(f)*as)));
% MDot     = n + ((b/(a*h*e))*(ar*((p*cos(f))-(2*r*e))-((p+r)*(sin(f*as)))));

orbDot = [aDot
          eDot
          iDot
          OmegaDot
          wDot
          fDot];

end