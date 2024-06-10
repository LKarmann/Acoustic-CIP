function [Anh, submesh] = inverse_solver_1(wh, mesh, h, s, Proposal)
% Solves the inverse problem using the knowledge of wh at internal nodes.
%
% Arguments:
% wh ('double'): Function which is the solution of Method 1.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% s ('scalar'): Pseudo-frequency at which equation is solved.
% Proposal ('int'): Proposal for a calculation at the vertices of the
%                   square.
%
% Returns:
% Anh ('double'): Reconstructed function on the internal mesh.
% submesh ('msh'): Mesh on which Anh has been calculated.
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
Wnh = wh(mesh.vtx(:,1) >=0 & mesh.vtx(:,2) >=0 & mesh.vtx(:,1) <=1 & mesh.vtx(:,2) <=1);

% Weighted mass matrix
Mnh = Wnh .* Mh;

% Boundary vector
Fnh = integral(Sigma2, Vh2, @(X) funFn(X, wh, mesh, submeshb, Proposal));

% Right-hand vector
bnh = (Fnh - Kh * Wnh) /(s^2);


% Solving the linear system
Anh = Mnh \ bnh;

end