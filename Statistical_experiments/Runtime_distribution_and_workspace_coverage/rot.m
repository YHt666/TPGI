function vr=rot(k,q,v)
% Rodrigues Rotation Formula.
% Input
% k：axis of rotation;
% q：angle of rotation;
% v：vector to be rotated.
% Output
% vr: vector obtained after rotation.
vr=cos(q)*v+(1-cos(q))*dot(k,v)*k+sin(q)*cross(k,v);