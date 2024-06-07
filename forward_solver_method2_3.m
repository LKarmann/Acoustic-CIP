function [U_list, mesh, t_list] = forward_solver_method2_3(h, dt, Nt, afun, pfun, om, delta, theta, Progression)
% Solves the wave equation with inhomogeneous Neumann-Absorbing boundary 
% conditions as explained in Method 2.3.
% See documentation of Gypsilab.
%
% Arguments:
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% dt ('scalar'): Time-step for the time discretization.
% Nt ('integer'): Number of time-steps.
% afun ('function_handle'): Wave speed function of 3D spatial variable.
% pfun ('function_handle'): Boundary condition function that appears in
%                           Method 2. Must be a function of time.
% delta ('scalar'): Parameter for the Newmark numerical integration.
%                   Must be between 1/2 and 1.
% theta ('scalar'): Parameter for the Newmark numerical integration.
%                   Must be between 0 and 1/2.
% Progression ('logical'): 1 means that the progression step is displayed.
%
% Returns:
% U_list (Nx(Nt+1) 'double'): Values of the solution u on the N nodes of
%                             the mesh at each time step.
% mesh ('msh'): Mesh used to solve.
%               See documentation of Gypsilab.
% t_list(1x(Nt+1) 'double'): Time discretization.



% Discretization of the time
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
bound_fun_lef = @(X) (abs(X(:,1)+0.5) < dx/3);
bound_fun_rig = @(X) (abs(X(:,1)-1.5) < dx/3);

B = integral(Sigma, Vh, bound_fun);

C1 = integral(Sigma, Vh, bound_fun_bot, Vh)+integral(Sigma, Vh, bound_fun_lef, Vh)+integral(Sigma, Vh, bound_fun_rig, Vh);
C2 = integral(Sigma, Vh, bound_fun, Vh);

C = @(t) C1 + (om * t > 2*pi)*C2;


U_list = damped_newmark(M,C,K,@(t) pfun(t)*B,zeros(size(B)),zeros(size(B)),dt,Nt,delta,theta,Progression);


end