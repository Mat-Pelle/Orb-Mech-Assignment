function orbDot = orbProp_gaussJ2(y, mu, j2, Earth_R)

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

%% Perturbations

%b = a*sqrt(1-(e^2));
p = a*(1-(e^2));%(b^2)/a;
r = p/(1+(e*cos(f)));
%v = sqrt(((2*mu)/r)-(mu/a));
%n = sqrt(mu/(a^3));
h = sqrt(p*mu);%n*a*b;

%--- J2 Perturbation                                                          
aj2_rsw = -(3/2)*((j2*mu*Earth_R^2)/r^4)*[1-(3*sin(i)^2*sin(f+w)^2)
                                             sin(i)^2*sin(2*(f+w))
                                             sin(2*i)*sin(f+w)];  
    
%Accelerations due to perturbations
ar = aj2_rsw(1);%+aMoon_rsw(1);
as = aj2_rsw(2);%+aMoon_rsw(2);
aw = aj2_rsw(3);%+aMoon_rsw(3);

%% Variation of keplerian elements (in RSW)

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