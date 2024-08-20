function y = refine_values(x,mesh,mesh_ref,method)
% Interpolates values of x defined on mesh to a refined mesh, mesh_ref.
% See documentation of Gypsilab.
%
% Arguments:
% x ('vector'): Values that must be interpolated.
% mesh ('msh'): Mesh on which x is given.
%               See documentation of Gypsilab.
% mesh_ref ('msh'): Refined mesh on which x must be interpolated.
%                   See documentation of Gypsilab.
% method ('str'): Method used for the interpolation. 
%                 'nearest': selects the first nearest point.
%                 'mean': makes an average of the values of the two nearest
%                 points.
%                 'linear': makes a linear approximation. In this method,
%                 it is assumed that mesh_ref is a refined mesh obtained by
%                 a midpoint algorithm.
%
% Returns:
% y ('vector'): Interpolated values of x.

y = zeros(size(mesh_ref.vtx(:,1)));

if method == "nearest"
    for k = 1:size(y,1)
        [~,I] = min((mesh.vtx(:,1)-mesh_ref.vtx(k,1)).^2 + (mesh.vtx(:,2)-mesh_ref.vtx(k,2)).^2);
    
        y(k) = x(I);
    end
elseif method == "mean"
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
elseif method == "linear"
    step2 = mesh_ref.stp;
    dx2 = step2(1);
    for k = 1:size(y,1)
        [M,I] = min((mesh.vtx(:,1)-mesh_ref.vtx(k,1)).^2 + (mesh.vtx(:,2)-mesh_ref.vtx(k,2)).^2);

        if M < dx2^2/5
            y(k) = x(I);
        else
            for j = 1:size(mesh.elt,1)
                d1 = ((mesh.vtx(mesh.elt(j,1),1) + mesh.vtx(mesh.elt(j,2),1))/2 - mesh_ref.vtx(k,1)).^2 ...
                    + ((mesh.vtx(mesh.elt(j,1),2) + mesh.vtx(mesh.elt(j,2),2))/2 - mesh_ref.vtx(k,2)).^2;
                d2 = ((mesh.vtx(mesh.elt(j,1),1) + mesh.vtx(mesh.elt(j,3),1))/2 - mesh_ref.vtx(k,1)).^2 ...
                    + ((mesh.vtx(mesh.elt(j,1),2) + mesh.vtx(mesh.elt(j,3),2))/2 - mesh_ref.vtx(k,2)).^2;
                d3 = ((mesh.vtx(mesh.elt(j,2),1) + mesh.vtx(mesh.elt(j,3),1))/2 - mesh_ref.vtx(k,1)).^2 ...
                    + ((mesh.vtx(mesh.elt(j,2),2) + mesh.vtx(mesh.elt(j,3),2))/2 - mesh_ref.vtx(k,2)).^2;
                if d1 < dx2^2/5
                    y(k) = (x(mesh.elt(j,1)) + x(mesh.elt(j,2)))/2;
                elseif d2 < dx2^2/5
                    y(k) = (x(mesh.elt(j,1)) + x(mesh.elt(j,3)))/2;
                elseif d3 < dx2^2/5
                    y(k) = (x(mesh.elt(j,2)) + x(mesh.elt(j,3)))/2;
                end
            end
        end
        if y(k) == 0
            disp('Error')
        end
    end
end
end