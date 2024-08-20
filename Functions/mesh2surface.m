function [x,y,z] = mesh2surface(u,mesh)
% Transforms a vector of the values associated to a mesh into a matrix form
% of a surface as usually in Matlab. mesh is assumed to be a standard mesh
% as given by functions like mshSquare or mshSquare2. Should not be used on
% an irregular mesh of a square.
%
% Arguments:
% u ((N^2)x1 'double'): Values of a function on the nodes of mesh. Indices
%                       between u and mesh should be coincident.
% mesh ('msh'): Mesh on which u is defined.
%               See documentation of Gypsilab.
%
% Returns:
% x (Nx1 'double'): x-discretisation coordinates of the square.
% y (Nx1 'double'): y-discretisation coordinates of the square.
% z (NxN 'double'): Values of u associated to the square mesh (x,y).


x = unique(mesh.vtx(:,1));
y = unique(mesh.vtx(:,2));

z = zeros(size(x,1));

for k = 1:size(z,2)
    z(:,k) = u(mesh.vtx(:,1)==x(k));
end

end