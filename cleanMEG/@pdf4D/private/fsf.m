function cost=fsf(xyzr,hs);

%FSF cost = fsf(xyzr,hs)
%    to be minimized in sphere fit
%
%    xyzr = [ctr;r]
%    where ctr - initial guess for center
%          r - initial guess for radius
%    hs - 3 by N head shape array

%1 by 3 sphere center
ctr = xyzr(1:3);
%sphere radius
r = xyzr(4);

%translate hs to move center to [0,0,0]
nhs = size(hs,2);
hs = hs - repmat(ctr,1,nhs);

%cost function
%sqrt(sum(hs.^2)) is 1 by nhs radii for all hs points
cost = sum(...
            (...
              sqrt(sum(hs.^2)) - r*ones(1,nhs) ...
             ).^2 ...
           );