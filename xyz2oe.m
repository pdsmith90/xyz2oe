function [a,e,i,w,OM,M] = xyz2oe(r,rdot)
global GM
% Question 3
% inverse transformation from state vector at time t to Keplerian oe

% angular momentum
h=cross(r,rdot);

% unit vectors
x=[1;0;0];
z=[0;0;1];

% n
n=cross(z,h);

% eccentricity vector
e3d=1./GM.*cross(rdot,h)-r./norm(r);

% eccentricity
e=norm(e3d);

% inclination
i=acos(dot(h,z)./(norm(h).*norm(z)));
if i>pi
    i=2.*pi-i;
end

% longitude of the ascending node
OM=acos(dot(x,n)./(norm(x).*norm(n)));
if n(2)<0
    OM-w.*pi-OM;
end

% argument of perigee
w=acos(dot(n,e3d)./(norm(n).*e));
if e3d(3)<0
    w=2.*pi-w;
end

% semi-major axis
a=norm(h).^2./(GM.*(1-e.^2));


% last element = many choices.

% true anomaly
f=acos(dot(e3d,r)./(e.*norm(r)));
if dot(r,rdot)<0
    f=2.*pi-f;
end

% cccentric anomaly
E=acos((e+cos(f))./(1+e.*cos(f)));

% mean anomaly
M=E-e.*sin(E);

end
