function mesh_ref = refine_mesh(mesh,vertices)
% Refines a mesh on every element having a vertex in vertices, using
% midpoint algorithm.
% See documentation of Gypsilab.
%
% Arguments:
% mesh ('msh'): Mesh which might be refined.
%               See documentation of Gypsilab.
% vertices ('vector'): List of the vertices in mesh which might be refined.
%
% Returns:
% mesh_ref ('msh'): Refined mesh.

l = false(size(mesh.elt(:,1)));

for k = 1:size(mesh.elt,1)

    l(k) = max(ismember(vertices,mesh.elt(k,:)));

end

list = find(l);

mesh_ref = mshMidpoint(mesh, list);

end