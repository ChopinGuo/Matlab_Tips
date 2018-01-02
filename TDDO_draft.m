syms sys1 sys2 u1 u2 x1 x2 x3 d1 d2 z1 z2
z_dot = [sys1; sys2];
z = [z1; z2];
x = [x1; x2; x3];
u = [u1; u2]

L = [40 0 -20; 0 -20 60];
A = [-0.8 1 1.6; 0 -3 2; 0 0 -6];
Bu = [0 0; 1 0; 0 1];
Bd = [0.8 0; 0 -1; -0.4 1.2];

z_dot = -L*Bd*(z+L*x)-L*(A*x+Bu*u)