function a_list = adaptive_gradient_1_1(u_til,a_0,gamma_0,q,tol,N,mesh_0,dt,Nt,pfun,method,bound,post_processing,intern)
% Reconstructs an approximation of a_* knowing only u on the boundary.
% It applies a simple gradient descent with adaptive mesh refinement.
% If a singular matrix appears in the solving, the values are set to an
% interpolated a_0.
% See Subsection 3.3.
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
% tol ('scalar'): Tolerance level for refining the mesh. Nodes with values
%                 greater than tol x max are refined. 
%                 Must be between 0 and 1.
% N ('int'): Number of iterations of the gradient method.
% mesh_0 ('msh'): Mesh used for the initialisation.
%                 See documentation of Gypsilab.
% dt ('scalar'): Time step for the time discretisation.
% Nt ('integer'): Number of time steps. Corresponds to Nt+1 time steps or
%                 T = Nt x dt.
% pfun ('function_handle'): Boundary condition function p(t). 
%                           Must be a function of time.
% method ('str'): Method used for the interpolation of the values at each
%                 refining. Should be "nearest" or "linear".
% bound ('logical'): Indicates if the scattered data are given on the whole
%                    boundary or only on the top side. 1 means that all the
%                    boundary is given.
% post_processing ('str'): Methods of post_processing the reconstructed
%                          values of a_n. Must be "None" or a method in:
%                          "Cutoff" means that values under 1 are set to 1.
%                          "Putup" means that values under 1 are set in the
%                          upright direction.
% intern ('logical'): 1 means that the values of the reconstructed a_n are
%                     set to 1 outside Omega_1.
%
% Returns:
% a_list (N+1 x 2 'cell'): Values of the reconstructed a_n at each
%                          iteration. Each row corresponds to an iteration.
%                          The first column contains the values of a_n and
%                          the second contains the mesh over which it is
%                          defined.


a_list = cell(N+1,2);
a_list{1,1} = a_0;
a_list{1,2} = mesh_0;


a_n = a_0;
mesh_n = mesh_0;
u_til_n = u_til;

step = mesh_0.stp;
dx = step(1);

if bound
    bound_mask_n = abs(mesh_0.vtx(:,1)+0.5) <= dx/3 | abs(mesh_0.vtx(:,2)+0.5) <= dx/3 | abs(mesh_0.vtx(:,1)-1.5) <= dx/3 | abs(mesh_0.vtx(:,2)-1.5) <= dx/3;
else
    bound_mask_n = abs(mesh_0.vtx(:,2)-1.5) <= dx/3;
end



for k = 1:N
    u_n = principal_solver_1_1(mesh_n,dt,Nt,a_n,pfun);

    deltau = u_til_n - u_n(bound_mask_n,:);

    lam_n = adjoint_solver_1_1(mesh_n,dt,Nt,a_n,deltau,bound);

    g_n = zeros(size(a_n));

    for j = 1:size(u_n,2)-1
        g_n = g_n - (u_n(:,j+1)-u_n(:,j)) .* (lam_n(:,j+1)-lam_n(:,j))/dt;
    end

    a_n = refine_values(a_0,mesh_0,mesh_n,"nearest") - k^q * g_n / gamma_0;

    if intern
        external = mesh_n.vtx(:,1) <= 0 | mesh_n.vtx(:,2) <= 0 | mesh_n.vtx(:,1) >= 1 | mesh_n.vtx(:,2) >= 1;
        a_n(external) = 1;
    end

    if post_processing == "Cutoff"
        a_n(a_n < 1) = 1;
    elseif post_processing == "Putup"
        a_n = 1 + abs(1-a_n);
    end

    if anynan(a_n)
        break
    end

    a_list{k+1,1} = a_n;
    a_list{k+1,2} = mesh_n;


    % Mesh refinement
    if k < N
    
        l = find(abs(1-a_n) >= tol*max(abs(1-a_n)));
    
        mesh_n_ref = refine_mesh(mesh_n, l);
    
        a_n = refine_values(a_n,mesh_n,mesh_n_ref,method);
    
        step = mesh_n_ref.stp;
        dx_ref = step(1);
        
        if bound
            bound_mask_n_ref = abs(mesh_n_ref.vtx(:,1)+0.5) <= dx_ref/3 | abs(mesh_n_ref.vtx(:,2)+0.5) <= dx_ref/3 | abs(mesh_n_ref.vtx(:,1)-1.5) <= dx_ref/3 | abs(mesh_n_ref.vtx(:,2)-1.5) <= dx_ref/3;
        else
            bound_mask_n_ref = abs(mesh_n_ref.vtx(:,2)-1.5) <= dx_ref/3;
        end
    
        u_til_n_ref = zeros([nnz(bound_mask_n_ref), size(u_til,2)]);
    
        for j = 1:size(u_til,2)
            u_til_n_ref(:,j) = refine_values_boundary(u_til_n(:,j),mesh_n.vtx(bound_mask_n,:),mesh_n_ref.vtx(bound_mask_n_ref,:),dx,method);
        end
    
        u_til_n = u_til_n_ref;
    
        dx = dx_ref;
        bound_mask_n = bound_mask_n_ref;
    
        mesh_n = mesh_n_ref;
    end
end

if anynan(a_n)
    a = zeros(size(mesh_n.vtx(:,1)));

    for j = k+1:N+1
        a_list{j,1} = a;
        a_list{j,2} = mesh_n;
    end
end
end