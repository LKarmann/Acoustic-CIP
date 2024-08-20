function y = refine_values_boundary(x,bound,bound_ref,dx,method)
% Interpolates values of x defined on the boundary of a mesh, bound, to the
% boundary of a refined mesh, bound_ref.
% See documentation of Gypsilab.
%
% Arguments:
% x ('vector'): Values that must be interpolated.
% bound ('double'): Coordinates of the nodes on the boundary of the mesh on
%                   which x is defined.
%                   See documentation of Gypsilab.
% bound_ref ('double'): Coordinates of the nodes on the boundary of the 
%                       refined mesh on which y is defined.
%                       See documentation of Gypsilab.
% dx ('scalar'): Minimal step size of the boundary bound.
% method ('str'): Method used for the interpolation. 
%                 'nearest' select the first nearest point.
%                 'linear' makes a linear approximation.
%
% Returns:
% y ('vector'): Interpolated values of x.

y = zeros(size(bound_ref(:,1)));

if method == "nearest"
    for k = 1:size(y,1)
        [~,I] = min((bound(:,1)-bound_ref(k,1)).^2 + (bound(:,2)-bound_ref(k,2)).^2);
    
        y(k) = x(I);
    end
elseif method == "linear"
    for k = 1:size(y,1)
        [M,I] = mink((bound(:,1)-bound_ref(k,1)).^2 + (bound(:,2)-bound_ref(k,2)).^2,2);

        if M(1) < dx^2/9
            y(k) = x(I(1));
        else
            if abs(bound(I(1),1)-bound(I(2),1)) >= abs(bound(I(1),2)-bound(I(2),2)) 
                y(k) = x(I(1)) + (x(I(2))-x(I(1)))*(bound_ref(k,1)-bound(I(1),1))/(bound(I(2),1)-bound(I(1),1));
            else
                y(k) = x(I(1)) + (x(I(2))-x(I(1)))*(bound_ref(k,2)-bound(I(1),2))/(bound(I(2),2)-bound(I(1),2));
            end
        end
    end
end
end