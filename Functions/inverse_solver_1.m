function [Ah,submesh] = inverse_solver_1(wh,mesh,h,s,Proposal)
% Solves Inverse problem 1 using the knowledge of wh at internal nodes.
% First implementation of the inverse solver, relies on various proposals
% for the computation of the normal derivative.
% See Appendix D.2.
%
% Arguments:
% wh ('double'): Laplace transform of the solution of a forward problem.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% s ('scalar'): Pseudo-frequency at which the equation is solved.
% Proposal ('int'): Proposal for a calculation at the vertices of the
%                   square.
%                   See Appendix D.2.
%
% Returns:
% Ah ('double'): Reconstructed acoustic coefficient.
% submesh ('msh'): Mesh of ]0, 1[x]0, 1[ on which Ah has been calculated.
%                  See documentation of Gypsilab.


N = ceil(2*sqrt(2)/h);                                % Number of points

while mod(N,4) ~= 0
    N = N+1;
end


% Creation of the submesh
submesh = mshSquare2(N/2, [0,1,0,1]);
submeshb = submesh.bnd;


% Integration domain
Omega2 = dom(submesh, 7);      % 1  3  7  12
Sigma2 = dom(submeshb, 4);     % 1  2  3  4  5


% Finite element
Vh2 = fem(submesh, 'P1');


% Stiffness matrix
Kh = integral(Omega2, grad(Vh2), grad(Vh2));

% Mass matrix
Mh = integral(Omega2, Vh2, Vh2);

% Weight vector
Wh = wh(mesh.vtx(:,1) >=0 & mesh.vtx(:,2) >=0 & mesh.vtx(:,1) <=1 & mesh.vtx(:,2) <=1);

% Weighted mass matrix
Mh_2 = Mh .* Wh';

% Boundary vector
Fh = integral(Sigma2, Vh2, @(X) funFn(X, wh, mesh, submeshb, Proposal));

% Right-hand vector
bh = (Fh - Kh * Wh) /(s^2);


% Solving the linear system
Ah = Mh_2 \ bh;

end