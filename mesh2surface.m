function [x,y,z] = mesh2surface(u,mesh)
% Transforms a vector of the values associated to a mesh into a matrix form
% of a surface as usually in Matlab. mesh is assumed to be a standard mesh
% as given by functions like mshSquare or mshSquare2. Should not be used on
% an irregular mesh of a square.
%
% Arguments:
% u (n^2x1 'double'): Values of a function on the nodes of mesh. Indices
%                     between u and mesh should be coincident.
% mesh ('msh'): Mesh on which u is defined.
%               See documentation of Gypsilab.
%
% Returns:
% x (nx1 'double'): x-discretization coordinates of the square.
% y (nx1 'double'): y-discretization coordinates of the square.
% z (nxn 'double'): Values of u associated to the square mesh (x,y).


x = unique(mesh.vtx(:,1));
y = unique(mesh.vtx(:,2));

z = zeros(size(x,1));

for k = 1:size(z,2)
    z(:,k) = u(mesh.vtx(:,1)==x(k));
end

end