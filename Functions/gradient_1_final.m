function a_list = gradient_1_final(u_til,a_0,gamma_0,q,tol,delta,N,mesh_0,dt,Nt,pfun,bound,mean,post_processing)
% Reconstructs an approximation of a_* knowing only u_til on the boundary.
% It applies a simple gradient descent. If a singular matrix appears in the
% solving, the values are set to a_0.
% Final version of the gradient method, values outside Omega_1 are set to 1.
% See Section 3 and Appendix J.
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
% tol ('scalar'): Tolerance parameter for the cut-off of the small values.
%                 Values lower than tol x max are set to 1.
%                 Should be 0 for no-post-processing.
% delta ('scalar'): Delta parameter for an improvement of the method. It
%                   consists in replacing the regularisation term by an
%                   integral over a small ring near the boundary of Omega
%                   with a width delta. Should be between 0 and 3/2.
%                   0 means that this method is not applied and we just
%                   apply the method described in Section 3.
%                   Not described in the report.
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
% mean ('logical'): Indicates if the mean post-processing is applied.
% post_processing ('str'): Methods of post_processing the reconstructed
%                          values of a_n. Must be "None" or a method in:
%                          "Smooth" means that the data are smoothed using
%                          a 2D-smoothing with a gaussian method.
%                          "Cutoff" means that values under 1 are set to 1.
%                          "Putup" means that values under 1 are set in the
%                          upright direction.
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

mask = abs(mesh_0.vtx(:,1)+0.5) <= delta | abs(mesh_0.vtx(:,2)+0.5) <= delta | abs(mesh_0.vtx(:,1)-1.5) <= delta | abs(mesh_0.vtx(:,2)-1.5) <= delta;


a_n = a_0;

for k = 1:N
    u_n = principal_solver_1_1(mesh_0,dt,Nt,a_n,pfun);

    deltau = u_til - u_n(bound_mask,:);

    lam_n = adjoint_solver_1_1(mesh_0,dt,Nt,a_n,deltau,bound);

    g_n = zeros(size(a_n));

    for j = 1:size(u_n,2)-1
        g_n = g_n - (u_n(:,j+1)-u_n(:,j)) .* (lam_n(:,j+1)-lam_n(:,j))/dt;
    end
    
    if delta > 0
        g_n(mask) = g_n(mask) + gamma_0 * (a_n(mask,:)-a_0(mask,:))/k^q;

        a_n = a_n - k^q * g_n / gamma_0;
    
    else
        a_n = a_0 - k^q * g_n / gamma_0;
    end

    if post_processing == "Cutoff"
        a_n(a_n < 1) = 1;
    elseif post_processing == "Putup"
        a_n = 1 + abs(1-a_n);
    end

    a_n(external) = 1;

    if anynan(a_n)
        break
    end

    if mean
        a_nbis = zeros(size(a_n));
        vertex_list = zeros(size(a_n));
    
        for j = 1:size(mesh_0.elt,1)
            x = (a_n(mesh_0.elt(j,1)) + a_n(mesh_0.elt(j,2)) + a_n(mesh_0.elt(j,3)))/3;
    
            vertex_list(mesh_0.elt(j,1)) = vertex_list(mesh_0.elt(j,1))+1;
            vertex_list(mesh_0.elt(j,2)) = vertex_list(mesh_0.elt(j,2))+1;
            vertex_list(mesh_0.elt(j,3)) = vertex_list(mesh_0.elt(j,3))+1;
            a_nbis(mesh_0.elt(j,1)) = a_nbis(mesh_0.elt(j,1)) + x;
            a_nbis(mesh_0.elt(j,2)) = a_nbis(mesh_0.elt(j,2)) + x;
            a_nbis(mesh_0.elt(j,3)) = a_nbis(mesh_0.elt(j,3)) + x;
        end

        a_n = a_nbis ./ vertex_list;
    end

    if tol > 0
        a_n(abs(a_n-1) < tol * max(abs(a_n-1))) = 1;
    end

    a_list(:,k+1) = a_n;
end

if anynan(a_n)
    for j = k+1:N+1
        a_list(:,j) = zeros(size(a_0));
    end
end
end