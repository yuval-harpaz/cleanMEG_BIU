function B=sarvas(r0,Q,r,n,ctr);
%SARVAS  B=sarvas(r0,Q,r,n,ctr)
%
%       Forward problem solution adapted from
%       Jukka Sarvas, Phys. Med. Biol., 1987, Vol 32, No 1, 11-22.
%
%   Input arguments:
%
%       r0  - 3 by nt  -  dipole position [m]
%                                    (nt is number of time points)
%       Q   - 3 by nt  -  dipole moment [Am]
%       r   - 3 by nch by nlp -  detector position [m]
%                                    (nch is number of detectors)
%                                    (nlp - number of coils per detector)
%       n   - 3 by nch by nlp -   normalized detector direction
%                                    (nch is number of detectors)
%                                    (nlp - number of coils per detector)
%       ctr - 3 by 1   -   center of sphere [m]
%
%   Output:
%
%       B   - nch by nt - magnetic field
%
%Copyleft 1996-2006, eugene.kronberg@uchsc.edu

mu0_over_4pi = 1.0e-7;		%mu0=4*pi*1e-7

%%%%%%%%%%%%%%%%%%%%%%%%%%% Argument Check %%%%%%%%%%%%%%%%%%
[three nch nlp]=size(r);
if three ~= 3
   disp('First dimension of ''r'' must be 3');
   return
end

[three nch_n nlp_n]=size(n);
if three ~= 3
   disp('First dimension of ''n'' must be 3');
   return
end

if nch ~= nch_n
   disp('Second dimensions of ''n'' and ''r'' must be equal');
   return
end

if nlp ~= nlp_n
   disp('Third dimensions of ''n'' and ''r'' must be equal');
   return
end

[three nt]=size(r0);
if three ~= 3
   disp('First dimension of ''r0'' must be 3');
   return
end

[three nt_q]=size(Q);
if three ~= 3
   disp('First dimension of ''Q'' must be 3');
   return
end

if nt ~= nt_q
   disp('Second dimension of ''r0'' and ''Q'' must be equal');
   return
end

[three one]=size(ctr);
if three ~= 3
   disp('First dimension of ''ctr'' must be 3');
   return
end
if one ~= 1
   disp('Second dimension of ''ctr'' must be 1');
   return
end

%%%%%%%%%%%%%%%%%%%%%%% Real Sarvas %%%%%%%%%%%%%%%%%%%%%

%all matrixes have dimensions 3 (x,y,z) by nch by nlp by nt:
r0=repmat(reshape(r0,3,1,1,nt),[1,nch,nlp,1]);
Q=repmat(reshape(Q,3,1,1,nt),[1,nch,nlp,1]);
r=repmat(r,[1,1,1,nt]);
n=repmat(n,[1,1,1,nt]);
ctr=repmat(ctr,[1,nch,nlp,nt]);

r0 = r0 - ctr;
r  = r  - ctr;
a  = r  - r0;

mod_a = sqrt(sum(a.*a));
mod_r = sqrt(sum(r.*r));
mod_r0 = sqrt(sum(r0.*r0));

ar = sum(a.*r);
r0r = sum(r0.*r);

F = repmat(mod_a.*(mod_r.*mod_a + mod_r.^2 - r0r),[3,1,1]);

gradF = ...
    repmat(mod_a.^2./mod_r + ar./mod_a + 2*mod_a + 2*mod_r,[3,1,1]).*r - ...
    repmat(mod_a + 2*mod_r + ar./mod_a,[3,1,1]).*r0;

Qxr0 = cross(Q, r0);
Qxr0r = repmat(sum(Qxr0.*r),[3,1,1]);

B = squeeze( ...
        sum( ...
           sum( ...
              (mu0_over_4pi*(F.*Qxr0 - Qxr0r.*gradF)./F.^2).*n ...
           ) ... %sum x,y,z
        ,3) ... %sum coils
    ); %final B should be nch by nt
