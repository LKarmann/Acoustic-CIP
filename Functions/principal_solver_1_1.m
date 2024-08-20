function u = principal_solver_1_1(mesh,dt,Nt,a,pfun)
% Solves the principal problem, Problem (1). Applies the method
% forward_solver_method2_1 with delta = 1/2 and theta = 1/2.
% See documentation of Gypsilab and Appendices C and I.
%
% Arguments:
% mesh ('msh'): Mesh used to solve the equation.
%               See documentation of Gypsilab.
% dt ('scalar'): Time step for the time discretisation.
% Nt ('integer'): Number of time steps. Corresponds to Nt+1 time steps or
%                 T = Nt x dt.
% a ('double'): Scattered values of an approximated acoustic wave
%               coefficient.
% pfun ('function_handle'): Boundary condition function p(t). 
%                           Must be a function of time.
%
% Returns:
% u (Nx(Nt+1) 'double'): Values of the solution uh on the N nodes of
%                        the mesh at each time step.


% Creation of the mesh
meshb = mesh.bnd;


% Integration domain
Omega = dom(mesh, 7);      % 1  3  7  12
Sigma = dom(meshb, 3);     % 1  2  3  4  5


% Finite element
Vh = fem(mesh, 'P1');


% Creation of auxiliar functions

function y = afun(X)

    y = zeros([size(X,1), 1]);

    for j = 1:size(y,1)
        [~,I] = min((mesh.vtx(:,1)-X(j,1)).^2 + (mesh.vtx(:,2)-X(j,2)).^2);
        y(j) = a(I);

    end
end



% Solving wave equation with Newmark method
K = integral(Omega, grad(Vh), grad(Vh));

M = integral(Omega, Vh, @afun, Vh);

step = mesh.stp;
dx = step(1);
bound_fun = @(X) (abs(X(:,2) - 1.5) < dx/3);

B = integral(Sigma, Vh, bound_fun);



% Undamped Newmark numerical integration

u = zeros([size(K,1), Nt+1]);

A = M + dt^2*K/2;

dA = decomposition(A);

for k = 3:Nt+1

    b = dt^2/2*(pfun((k-1)*dt) + pfun((k-3)*dt))*B + 2*M*u(:,k-1);


    u(:,k) = dA\b - u(:,k-2);

end
end