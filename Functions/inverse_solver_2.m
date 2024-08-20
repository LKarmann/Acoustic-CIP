function [Ah,submesh] = inverse_solver_2(wh,mesh,h,s,pfun,Extraction)
% Solves Inverse problem 1 using the knowledge of wh at internal nodes.
% Uses the method of reconstruction in Omega instead of Omega_1.
% See Appendix D.3.
%
% Arguments:
% wh ('double'): Laplace transform of the solution of a forward problem.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% s ('scalar'): Pseudo-frequency at which the equation is solved.
% pfun ('function_handle'): Boundary condition function that appears in
%                           Method 1. Should be a function of the space.
%                           Corresponds to the Laplace transform of p(t).
% Extraction ('logical'): 1 means that Ah are the extracted values in
%                         Omega_1. 0 means that Ah are the values in Omega.
%
% Returns:
% Ah ('double'): Reconstructed acoustic coefficient.
% submesh ('msh'): Mesh on which Ah has been calculated.
%                  If Extraction = 1, then submesh is a mesh of 
%                  ]0, 1[x]0, 1[. Else, submesh is mesh.
%                  See documentation of Gypsilab.


% Creation of the boundary
meshb = mesh.bnd;


% Integration domain
Omega = dom(mesh, 7);      % 1  3  7  12
Sigma = dom(meshb, 3);     % 1  2  3  4  5


% Finite element
Vh = fem(mesh, 'P1');


% Stiffness matrix
Kh = integral(Omega, grad(Vh), grad(Vh));

% Mass matrix
Mh = integral(Omega, Vh, Vh);

% Weight vector
Wh = wh;

% Weighted mass matrix
Mh_2 = Mh .* Wh';

% Boundary vector
Fh = integral(Sigma, Vh, pfun);

% Right-hand vector
bh = (Fh - Kh * Wh) /(s^2);


% Solving the linear system
Ah = Mh_2 \ bh;

if Extraction

    N = sqrt(mesh.size(1)/2);                    % Parameter of the mesh
    
    submesh = mshSquare2(N/2, [0,1,0,1]);

    Ah = Ah(mesh.vtx(:,1) >=0 & mesh.vtx(:,2) >=0 & mesh.vtx(:,1) <=1 & mesh.vtx(:,2) <=1);

else
    submesh = mesh;
end

end