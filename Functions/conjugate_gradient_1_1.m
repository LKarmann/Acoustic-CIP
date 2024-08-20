function a_list = conjugate_gradient_1_1(u_til,a_0,gamma_0,q,N,mesh_0,dt,Nt,pfun,bound,post_processing,intern)
% Reconstructs an approximation of a_* knowing only u_til on the boundary.
% It applies a conjugate gradient descent. If a singular matrix appears in
% the solving, the values are set to a_0.
% See Appendix J.2.
%
% Arguments:
% u_til (Mx(Nt+1) 'double'): Scattered values of the measurement on the
%                            boundary. Could be on the whole boundary or
%                            only on the top side.
% a_0 ('double'): Scattered values of the initial guess for the wave
%                 coefficient.
% gamma_0 ('scalar'): Regularisation parameter for the Tikhonov functional.
% q ('scalar'): Exponent parameter for the decreasing of the regularisation
%               parameter. Should be between 0 and 1.
% N ('int'): Number of iterations of the gradient method.
% mesh_0 ('msh'): Mesh used to solve the equation.
%                 See documentation of Gypsilab.
% dt ('scalar'): Time step for the time discretisation.
% Nt ('integer'): Number of time steps. Corresponds to Nt+1 time steps or
%                 T = Nt x dt.
% pfun ('function_handle'): Boundary condition function p(t). 
%                           Must be a function of time.
% bound ('logical'): Indicates if the scattered data are given on the whole
%                    boundary or only on the top side. 1 means that all the
%                    boundary is given.
% post_processing ('str'): Methods of post_processing the reconstructed
%                          values of a_n. Must be "None" or a method in:
%                          "Smooth" means that the data are smoothed using
%                          a 2D-smoothing with a gaussian method.
%                          "Cutoff" means that values under 1 are set to 1.
%                          "Putup" means that values under 1 are set in the
%                          upright direction.
% intern ('logical'): 1 means that the values of the reconstructed a_n are
%                     set to 1 outside Omega_1.
%
% Returns:
% a_list ('double'): Values of the reconstructed a_n at each
%                    iteration on mesh.


a_list = zeros(size(a_0,1),N+1);
a_list(:,1) = a_0;

step = mesh_0.stp;
dx = step(1);

if bound
    bound_mask = abs(mesh_0.vtx(:,1)+0.5) <= dx/3 | abs(mesh_0.vtx(:,2)+0.5) <= dx/3 | abs(mesh_0.vtx(:,1)-1.5) <= dx/3 | abs(mesh_0.vtx(:,2)-1.5) <= dx/3;
else
    bound_mask = abs(mesh_0.vtx(:,2)-1.5) <= dx/3;
end

external = mesh_0.vtx(:,1) <= 0 | mesh_0.vtx(:,2) <= 0 | mesh_0.vtx(:,1) >= 1 | mesh_0.vtx(:,2) >= 1;


a_n = a_0;

norm_g_nm1 = 1;
d_n = zeros(size(a_0));

for k = 1:N
    u_n = principal_solver_1_1(mesh_0,dt,Nt,a_n,pfun);

    deltau = u_til - u_n(bound_mask,:);

    lam_n = adjoint_solver_1_1(mesh_0,dt,Nt,a_n,deltau,bound);

    g_n = zeros(size(a_n));

    for j = 1:size(u_n,2)-1
        g_n = g_n - (u_n(:,j+1)-u_n(:,j)) .* (lam_n(:,j+1)-lam_n(:,j))/dt;
    end

    g_n = g_n + gamma_0/k^q * (a_n - a_0);

    beta = norm(g_n)^2 / norm_g_nm1^2;

    d_n = -g_n + beta * d_n;

    alpha = (g_n.' * d_n)/norm(d_n)^2;

    a_n = a_0 - k^q * alpha * d_n / gamma_0;

    norm_g_nm1 = norm(g_n);

    if post_processing == "Cutoff"
        a_n(a_n < 1) = 1;
    elseif post_processing == "Putup"
        a_n = 1 + abs(1-a_n);
    elseif post_processing == "Smooth"
        [~, ~, b_n] = mesh2surface(a_n,mesh_0);
        b_n = smoothdata2(b_n,'gaussian');
        a_n = reshape(b_n,[],1);
    end

    if intern
        a_n(external) = 1;
    end

    if anynan(a_n)
        break
    end

    a_list(:,k+1) = a_n;
end

if anynan(a_n)
    for j = k+1:N+1
        a_list(:,j) = zeros(size(a_0));
    end
end
end