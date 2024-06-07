function A = boundary_post_processing(Anh, submesh)
% Realises a post-processing on the boundary of the reconstructed values of
% Anh. Consists in setting to 1 the boundary values of Anh. Anh is supposed
% to be on a submesh of ]0;1[x]0;1[.
%
% Arguments:
% Anh ('double'): The values of the reconstructed Anh over the submesh.
% submesh ('msh'): The mesh of the internal square over which Anh has been
%                  calculated.
%                  See documentation of Gypsilab.
%
% Returns:
% A ('double'): The values of Anh post-processed. The values are unchanged
%               inside the square and set to 1 on the boundary.

step = submesh.stp;
dx = step(1);

B = (submesh.vtx(:,1) >= dx/3 & submesh.vtx(:,1) <= (1-dx/3) & submesh.vtx(:,2) >= dx/3 & submesh.vtx(:,2) <= (1-dx/3));

A = Anh .* B + (1-B);

end