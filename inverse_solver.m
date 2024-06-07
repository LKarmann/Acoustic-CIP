function [Anh, submesh] = inverse_solver(wh, mesh, s)
% Solves the inverse problem using the knowledge of wh at internal nodes.
%
% Arguments:
% wh ('double'): Function which is the solution of Method 1.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% s ('scalar'): Pseudo-frequency at which equation is solved.
%
% Returns:
% Anh ('double'): Reconstructed function on the internal mesh.
% submesh ('msh'): Mesh on which Anh has been calculated.
%                  See documentation of Gypsilab.


N = round(sqrt(size(mesh.vtx,1))-1);                                % Number of points


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

% Auxiliar function
function y = auxfunFn(X, wh, mesh, submeshb)
% Calculates the normal derivative of wh on the boundary using discrete
% derivative.
%
% Arguments:
% X (nx3 'double'): Position for a calculation of an integral.
%                   See documentation of Gypsilab.
% wh ('double'): Function which is the solution of Method 1.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% submeshb ('msh'): Mesh of the boundary of the sub-open set.
%                   See documentation of Gypsilab.
%
% Returns:
% y ('double'): Value which is used for the integration.
%               See documentation of Gypsilab.


n = sqrt(mesh.size(1)/2);                    % Parameter of the mesh

fn = zeros([size(submeshb,1) 1]);

dx = mesh.vtx(n+2,1) - mesh.vtx(1,1);
dy = mesh.vtx(2,2) - mesh.vtx(1,2);

for k = 1:size(submeshb,1)
    x = submeshb.vtx(k,1);
    y = submeshb.vtx(k,2);

    [~,l] = min((mesh.vtx(:,1)-x).^2 + (mesh.vtx(:,2)-y).^2);

    if x == 0 && y > 0 && y < 1
        fn(k) = (wh(l-n-1) - wh(l+n+1))/(2*dx);
    elseif x == 1 && y > 0 && y < 1
        fn(k) = (wh(l+n+1) - wh(l-n-1))/(2*dx);
    elseif y == 0 && x > 0 && x < 1
        fn(k) = (wh(l-1) - wh(l+1))/(2*dy);
    elseif y == 1 && x > 0 && x < 1
        fn(k) = (wh(l+1) - wh(l-1))/(2*dy);
    end

end


v = submeshb.vtx;


% Proposal 3: Use the nearest node
v(v(:,1) == 0 & v(:,2) == 0, 1:2) = [2,2];
v(v(:,1) == 0 & v(:,2) == 1, 1:2) = [2,2];
v(v(:,1) == 1 & v(:,2) == 0, 1:2) = [2,2];
v(v(:,1) == 1 & v(:,2) == 1, 1:2) = [2,2];


y = zeros([size(X,1) 1]);

for k = 1:size(y,1)
    [~, I] = min((v(:,1)-X(k,1)).^2 + (v(:,2)-X(k,2)).^2);
    y(k) = fn(I);
end

end


Fnh = integral(Sigma2, Vh2, @(X) auxfunFn(X, wh, mesh, submeshb));

% Right-hand vector
bnh = (Fnh - Kh * Wnh) /(s^2);


% Solving the linear system
Anh = Mnh \ bnh;

end