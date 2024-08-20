function [Uh,mesh,t_list] = forward_solver_method2_2(h,dt,Nt,afun,pfun,t_1,delta,theta)
% Solves the wave equation with mixed Neumann-Absorbing boundary conditions 
% as described in Problem (2). Uses the damped Newmark numerical 
% integration to solve the problem.
% See documentation of Gypsilab and Appendix C.
%
% Arguments:
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% dt ('scalar'): Time step for the time discretisation.
% Nt ('integer'): Number of time steps. Corresponds to Nt+1 time steps or
%                 T = Nt x dt.
% afun ('function_handle'): Wave speed function of 3D spatial variable.
% pfun ('function_handle'): Boundary condition function p(t). 
%                           Must be a function of time.
% t_1 ('scalar'): Time instant for changing the boundary condition on C(t).
% delta ('scalar'): Parameter for the Newmark numerical integration.
%                   Must be between 1/2 and 1.
% theta ('scalar'): Parameter for the Newmark numerical integration.
%                   Must be between 0 and 1/2.
%
% Returns:
% Uh (Nx(Nt+1) 'double'): Values of the solution uh on the N nodes of
%                         the mesh at each time step.
% mesh ('msh'): Mesh used to solve.
%               See documentation of Gypsilab.
% t_list (1x(Nt+1) 'double'): Time discretisation.


% Discretisation of the time
t_list = 0:dt:Nt*dt;


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


% Solving wave equation with Newmark method
K = integral(Omega, grad(Vh), grad(Vh));

M = integral(Omega, Vh, afun, Vh);

step = mesh.stp;
dx = step(1);
bound_fun = @(X) (abs(X(:,2) - 1.5) < dx/3);
bound_fun_bot = @(X) (abs(X(:,2)+0.5) < dx/3);

B = integral(Sigma, Vh, bound_fun);

C1 = integral(Sigma, Vh, bound_fun_bot, Vh);
C2 = integral(Sigma, Vh, bound_fun, Vh);


Uh = damped_newmark(M,C1,C2,K,@(t) pfun(t)*B,zeros(size(B)),zeros(size(B)),dt,Nt,t_1,delta,theta);

end