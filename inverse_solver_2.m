function [Anh, submesh] = inverse_solver_2(wh, mesh, h, s, pfun, Extraction)
% Solves the inverse problem using the knowledge of wh at internal nodes
% according to Idea II.2.
%
% Arguments:
% wh ('double'): Function which is the solution of Method 1.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% s ('scalar'): Pseudo-frequency at which equation is solved.
% pfun ('function_handle'): Boundary condition function that appears in
%                           Method 1.
% Extraction ('logical'): 1 means that the return is the extracted value of
%                         Anh over the internal subset.
%
% Returns:
% Anh ('double'): Reconstructed function on the mesh.
% submesh ('msh'): Mesh on which Anh has been calculated.
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
Wnh = wh;

% Weighted mass matrix
Mnh = Wnh .* Mh;

% Boundary vector
Fh = integral(Sigma, Vh, pfun);

% Right-hand vector
bnh = (Fh - Kh * Wnh) /(s^2);


% Solving the linear system
Anh = Mnh \ bnh;

if Extraction

    N = sqrt(mesh.size(1)/2);                    % Parameter of the mesh
    
    submesh = mshSquare2(N/2, [0,1,0,1]);

    Anh = Anh(mesh.vtx(:,1) >=0 & mesh.vtx(:,2) >=0 & mesh.vtx(:,1) <=1 & mesh.vtx(:,2) <=1);

else
    submesh = mesh;
end

end