function A = boundary_post_processing(Ah,submesh)
% Realises a post-processing on the boundary of the reconstructed values of
% Ah. Consists in setting to 1 the boundary values of Ah. Ah is supposed to
% be defined on a mesh of ]0, 1[x]0, 1[, submesh.
%
% Arguments:
% Ah ('double'): Values of the reconstructed Ah over the mesh.
% submesh ('msh'): Mesh of ]0, 1[x]0, 1[ over which Ah has been calculated.
%                  See documentation of Gypsilab.
%
% Returns:
% A ('double'): Values of Ah after post-processing. Values are unchanged
%               inside the square and set to 1 on the boundary.

step = submesh.stp;
dx = step(1);

B = (submesh.vtx(:,1) >= dx/3 & submesh.vtx(:,1) <= (1-dx/3) & submesh.vtx(:,2) >= dx/3 & submesh.vtx(:,2) <= (1-dx/3));

A = Ah .* B + (1-B);

end