function [Anh, submesh, s, error_k] = blind_inverse_solver(wh, mesh, h, s_inf, s_sup, s_step, tol, pfun, Order)
% Solves the inverse problem without the knowledge of the pseudo-frequency
% s_star. Gives the best approximation of the solution Anh according to the
% relative quadratic error with wh.
%
% Arguments:
% wh ('double'): Function which is the solution of Method 1 for a value 
%                s_star.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% h ('scalar'): Mesh size parameter. Mesh grid should not have a finite
%               element diameter greater than h.
% s_inf ('scalar'): Inf value of s.
% s_sup ('scalar'): Sup value of s.
% s_step ('scalar'): Step for the discretization in s.
% tol ('scalar'): Tolerance to stop the loop.
% pfun ('function_handle'): Boundary condition function that appears in
%                           Method 1.
% Order ('string'): Order of visiting of the values of s. "increase" for
%                   an increasing visit, "decrease" otherwise.
%
% Returns:
% Anh ('double'): Reconstructed function on the internal mesh.
% submesh ('msh'): Mesh on which Anh has been calculated.
%                  See documentation of Gypsilab.
% s ('scalar'): Approximated value of the value of s_star.
% error_k ('scalar'): Relative quadratic error on wh for the reconstructed
%                     function Anh.


if Order == "increase"
    s_list = s_inf:s_step:s_sup;
else
    s_list = s_sup:-s_step:s_inf;
end


% Auxiliar function
function Y = Anfun(X, Anh, submesh)
% Auxiliar function for solving forward problem with the approximation of a
% given by Anh over submesh.

Y = zeros([size(X,1), 1]);

for j = 1:size(X,1)
    if min(X(j,1), X(j,2)) > 0 && max(X(j,1), X(j,2)) < 1

        [~,I] = min((submesh.vtx(:,1)-X(j,1)).^2 + (submesh.vtx(:,2)-X(j,2)).^2);

        Y(j) = Anh(I);

    else

        Y(j) = 1;

    end
end

end



% Initialising

s = s_list(1);

[Anh, submesh] = inverse_solver(wh, mesh, s);


% Calculation of the relative error on wh

wnh = forward_solver_method1(h, s, @(X) Anfun(X,Anh,submesh), pfun);

error_k = sqrt(sum((wnh-wh).^2))/sqrt(sum(wh.^2));

error_km1 = 2*error_k;



k = 2;

while k <= size(s_list,2) && error_k <= error_km1 && error_k > tol

    s_km1 = s;
    Anm1h = Anh;

    s = s_list(k);

    [Anh, submesh] = inverse_solver(wh, mesh, s);


    % Calculation of the relative error on wh

    wnh = forward_solver_method1(h, s, @(X) Anfun(X,Anh,submesh), pfun);

    error_km1 = error_k;

    error_k = sqrt(sum((wnh-wh).^2))/sqrt(sum(wh.^2));

    k = k+1;

end



if error_k > error_km1
    s = s_km1;
    Anh = Anm1h;
    error_k = error_km1;
end
end