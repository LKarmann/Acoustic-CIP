function [wh, mesh] = forward_solver_method1(h, s, afun, pfun)
% Solves the screened Poisson equation as explained in Method 1.
% See documentation of Gypsilab.
%
% Arguments:
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% s ('scalar'): Pseudo-frequency at which equation is solved.
% afun ('function_handle'): Wave speed function of 3D spatial variable.
% pfun ('function_handle'): Boundary condition function that appears in
%                           Method 1.
%
% Returns:
% wh (Nx1 'double'): Values of the solution w on the N nodes of the mesh.
% mesh ('msh'): Mesh used to solve.
%               See documentation of Gypsilab.


% Creation of the mesh
N = ceil(2*sqrt(2)/h);                                % Number of points, should be divisible by 4.

while mod(N,4) ~= 0
    N = N+1;
end

mesh = mshSquare2(N, [-0.5 1.5 -0.5 1.5]);
meshb = mesh.bnd;


% Integration domain
Omega = dom(mesh, 7);      % 1  3  7  12
Sigma = dom(meshb, 3);     % 1  2  3  4  5


% Finite element
Vh = fem(mesh, 'P1');


% Solving screened Poisson problem
K = integral(Omega, grad(Vh), grad(Vh)) + s^2 * integral(Omega, Vh, afun, Vh);

F = integral(Sigma, Vh, pfun);

wh = K \ F;

end