function M_rotate = quadRotation(phi)

%Computes a 4x4 matrix for rotating quadratures of A and B, to be
%multiplied by Theta functions, for V, need to square everything.

M_rotate = [cos(phi) sin(phi) 0 0;-sin(phi) cos(phi) 0 0;...
    0 0 cos(2.*phi) sin(2.*phi);0 0 -sin(2.*phi) cos(2.*phi)];
