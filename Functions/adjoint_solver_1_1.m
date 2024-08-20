function lam = adjoint_solver_1_1(mesh,dt,Nt,a,deltau,bound)
% Solves the adjoint problem, Problem (10). Applies a similar method to
% forward_solver_method2_1 with delta = 1/2 and theta = 1/2.
% See documentation of Gypsilab and Appendix I.
%
% Arguments:
% mesh ('msh'): Mesh used to solve the equation.
%               See documentation of Gypsilab.
% dt ('scalar'): Time step for the time discretisation.
% Nt ('integer'): Number of time steps. Corresponds to Nt+1 time steps or
%                 T = Nt x dt.
% a ('double'): Scattered values of an approximated acoustic wave
%               coefficient.
% deltau ('double'): Scattered values of u_til - u on a part of the
%                    boundary. Might be either the whole boundary or the
%                    top side of the square.
% bound ('logical'): Indicates if the scattered data are given on the whole
%                    boundary or only on the top side. 1 means that all the
%                    boundary is given.
%
% Returns:
% lam (Nx(Nt+1) 'double'): Values of the solution lambda on the N nodes of
%                          the mesh at each time step.


% Creation of the mesh
meshb = mesh.bnd;


% Integration domain
Omega = dom(mesh, 7);      % 1  3  7  12
Sigma = dom(meshb, 3);     % 1  2  3  4  5


% Finite element
Vh = fem(mesh, 'P1');


% Creation of auxiliar functions

step = mesh.stp;
dx = step(1);

function y = afun(X)

    y = zeros([size(X,1), 1]);

    for j = 1:size(y,1)
        [~,I] = min((mesh.vtx(:,1)-X(j,1)).^2 + (mesh.vtx(:,2)-X(j,2)).^2);
        y(j) = a(I);

    end
end

boundaryfun = @(X) (abs(X(:,2) - 1.5) < dx/3);


% Solving wave equation with Newmark method
K = integral(Omega, grad(Vh), grad(Vh));

M = integral(Omega, Vh, @afun, Vh);

if bound
    B = integral(Sigma, Vh, Vh);
else
    B = integral(Sigma, Vh, boundaryfun, Vh);
end

Deltau = zeros([size(B,1), size(deltau,2)]);

if bound
    Deltau(abs(mesh.vtx(:,1)+0.5) <= dx/3 | abs(mesh.vtx(:,2)+0.5) <= dx/3 | abs(mesh.vtx(:,1)-1.5) <= dx/3 | abs(mesh.vtx(:,2)-1.5) <= dx/3,:) = deltau;
else
    Deltau(abs(mesh.vtx(:,2)-1.5) <= dx/3,:) = deltau;
end


% Undamped Newmark numerical integration

lam = zeros([size(K,1), Nt+1]);

A = M + dt^2*K/2;

dA = decomposition(A);


for k = Nt-1:-1:1

    b = dt^2*B*(Deltau(:,k) + Deltau(:,k+2))/2 + 2*M*lam(:,k+1);


    lam(:,k) = dA\b - lam(:,k+2);

end
end