
function kep = UT_car2kep(r,v,mu)

rmod = norm(r);
vmod = norm(v);
vr = dot(r,v);

%total energy calculation
%xi = ((vmod^2)/2)-(mu/rmod);

%semi-mayor axis calculation
%a = -(mu)/(2*xi);
a = -mu / (vmod^2 - 2*mu/rmod);

%Angular momentum vector
h = cross(r,v);
hmod = norm(h);

%Inclination i
i = acos(h(3)/hmod);
  
%Ascending node n
k = [0 0 1];
n = cross(k,h)/norm(cross(h,k));
nmod = norm(n);

%Longitude of ascending node
if n(2)>=0
    omega = acos(n(1));
else
    omega = (2*pi)-acos(n(1));
end

%Excentricity e 
e = (1/mu)*cross(v,h)-r/rmod;
emod = norm(e); 

%Argument of periapsis
if e(3)>= 0
    w = acos(dot(n,e)/(nmod*emod));
else
    w = (2*pi)-acos(dot(n,e)/(nmod*emod));
end

%True anomaly
if vr>=0
    f = acos(dot(e,r)/(emod*rmod));
%     f = unwrap(f);
%     f = rad2deg(f);
else
    f = (2*pi) - acos(dot(e,r)/(emod*rmod)); 
%     f = unwrap(f);
%     f = rad2deg(f);
end

kep = [a emod i omega w f];
end
    

    