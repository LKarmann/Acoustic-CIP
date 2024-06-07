function y = refine_values(x,mesh,mesh_ref,method)
% Interpolates values of x on mesh to a refined mesh, mesh_ref.
% See documentation of Gypsilab.
%
% Arguments:
% x ('vector'): Values that must be interpolated.
% mesh ('msh'): Mesh on which x is given.
%               See documentation of Gypsilab.
% mesh_ref ('msh'): Refined mesh on which x must be interpolated.
%                   See documentation of Gypsilab.
% method ('str'): Method used for the interpolation. 
%                 'nearest' select the first nearest point.
%                 'linear' makes a linear approximation. In this method, it
%                 is assumed that mesh_ref is a refined mesh obtained by a
%                 midpoint algorithm.
%
% Returns:
% y ('vector'): Interpolated values of x.

y = zeros(size(mesh_ref.vtx(:,1)));

if method == "nearest"
    for k = 1:size(y,1)
        [~,I] = min((mesh.vtx(:,1)-mesh_ref.vtx(k,1)).^2 + (mesh.vtx(:,2)-mesh_ref.vtx(k,2)).^2);
    
        y(k) = x(I);
    end
elseif method == "linear"
    step = mesh.stp;
    dx = step(1);
    for k = 1:size(y,1)
        [M,I] = mink((mesh.vtx(:,1)-mesh_ref.vtx(k,1)).^2 + (mesh.vtx(:,2)-mesh_ref.vtx(k,2)).^2,4);

        if M(1) < dx^2/5
            y(k) = x(I(1));
        else
            y(k) = sum(x(I(M < dx^2 * 0.6)))/sum(M < dx^2);
        end
    end
end
end