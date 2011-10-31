function [ctr,r]=s_fit(pdf);

%S_FIT [ctr,r]=s_fit(pdf)
%    Sphere fit
%
%    pdf - pdf4D object
%    ctr - 3 by 1 - center of the sphere [cm]
%    r - radius of the sphere [cm]

%get headshape from pdf [m]
hs = get(pdf, 'hspoints');

%min and max for x (fron and back of the head)
xmin=min(hs(1,:));
xmax=max(hs(1,:));
%radius estimation for x
xradius=(xmax-xmin)/2;

%min and max for y (left and right of the head)
ymin=min(hs(2,:));
ymax=max(hs(2,:));
%radius estimation for y
yradius=(ymax-ymin)/2;

%initial value for the sphere radius
radius=mean([xradius yradius]);

%max for z (top of the head)
zmax=max(hs(3,:));

%initial value for the sphere center
ctr=[mean([xmin xmax]); ...
     mean([ymin ymax]); ...
     zmax-radius];

%Nelder-Mead (aka simplex)
xyzr = fminsearch(@(x) fsf(x,hs), [ctr;radius]);

%convert result to cm
ctr = 100*xyzr(1:3);
r = 100*xyzr(4);