function [pos, dir] = pos2arr(pdf, pos_in)

%POS2ARR   [pos, dir] = pos2arr(pdf, pos_in)
%Convert from position structure to pos and dir arrays
%
%pos_in is array structures with fields "position" and "direction"
%pos is 3 by nch by ncl array of coil position
%dir is 3 by nch by ncl array of coil direction
%nch is number of channels
%ncl is number coils per channel

nch = numel(pos_in);
%all channels should be of the same kind (same number of coils)!!!
ncl = size(pos_in(1).position,2);

%init arrays
pos = zeros(3,nch,ncl);
dir = zeros(3,nch,ncl);

for ch = 1:nch
    pos(:,ch,:) = pos_in(ch).position;
    dir(:,ch,:) = pos_in(ch).direction;
end
    