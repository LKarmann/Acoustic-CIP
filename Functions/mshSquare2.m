function mesh = mshSquare2(N,L)
% Creates a standard mesh on a square using Delaunay Triangulation.
% Adaptation of the function mshSquare from Gypsilab.
%
% Arguments:
% N ('int'): Number of triangles by side (N+1 points by side).
% L (4x1 'double'): Position of the vertices of the square.
%                   If L = [a b c d], then the square is [a, b]x[c, d].
%
% Returns:
% mesh ('msh'): Mesh of the class 'msh' from "openMsh" in Gypsilab.
%               See documentation of Gypsilab.


len = [L(2)-L(1), L(4)-L(1)];

% Delaunay mesh
x = L(1) : len(1)/N : L(2);
y = L(3) : len(2)/N : L(4);

% Triangulation
[x,y]  = meshgrid(x,y);
z      = zeros(size(x));
DT     = delaunayTriangulation([x(:) y(:)]);
Points = [DT.Points z(:)];

% Build mesh
mesh = msh(Points,DT.ConnectivityList);
end