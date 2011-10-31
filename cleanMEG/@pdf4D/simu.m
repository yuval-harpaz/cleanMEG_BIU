function [b, chl] = simu(obj, r, q, ctr)

%SIMU Data simulator [b, chl] = simu(obj, r, q, ctr)
%obj - pdf4D object
%r   - 3 by nt (number of time points) array of dipole position [cm]
%q   - 3 by nt array of dipole strength [nAm]
%ctr - 3 by 1 vector of sphere center [cm]
%
%Outputs:
%b   - nch (number of channels) by nt magnetic field
%chl - cell array of labels for meg channels

%test if pdf still exist
test(obj)

%check number of inputs
if nargin ~= 4
    error('Wrong number of inputs');
end

%get meg channel inexes
chi = channel_index(obj, 'meg', 'name');

%no meg channels - no simulation
if isempty(chi)
    return
end

%meg channel labels
chl = channel_label(obj, chi);

%get meg channel position and direction
ch = channel_position(obj, chi);

%convert position and direction for the "sarvas"
for c = 1:length(ch)
    pos(:,c,:) = ch(c).position;
    dir(:,c,:) = ch(c).direction;
end

%forward solution (scale to m and Am)
%(calls function in pdf4D/private)
b = sarvas(r*10^-2, q*10^-9, pos, dir, ctr*10^-2);