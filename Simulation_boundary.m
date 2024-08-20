%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Simulations and figures of the results of the reconstruction using
% boundary measurements and gradient method.


%% Preamble

clear;
close all;
format long;

addpath \gypsilab-master\gypsilab-master\openMsh;
addpath \gypsilab-master\gypsilab-master\openDom;
addpath \gypsilab-master\gypsilab-master\openFem;
addpath \Functions;


%% Definition of the parameters

h = 2^(-5);                                             % Mesh size

dt = 0.001;                                             % Time step (in seconds)

Nt = 2000;                                              % Number of time steps

t_list = 0:dt:Nt*dt;                                    % Time discretization

om = 80;                                                % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source

% Test 1: 1 inclusion
afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% Test 2: 2 inclusions
% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

delta = 1/2;                                            % Parameter of the Newmark method

theta = 1/2;                                            % Parameter of the Newmark method

N = 20;                                                 % Number of iterations for the gradient descent

gamma_0 = 0.01;                                         % Regularisation parameter

gamma_0bis = 1;                                         % Regularisation parameter for no post-processing

q = 1;                                                  % Regularisation exponent

qbis = 0.1;                                             % Regularisation exponent for no post-processing

tol = 0.95;                                             % Tolerance parameter

method = "linear";                                      % Interpolation for the refinement


%% Generation of the measurements

[u_list1, mesh_0, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);

step = mesh_0.stp;
dx = step(1);

% Extraction of the top side of the domain
bound_mask = abs(mesh_0.vtx(:,2)-1.5) <= dx/3;

% Extraction of the whole boundary
% bound_mask = abs(mesh_0.vtx(:,1)+0.5) <= dx/3 | abs(mesh_0.vtx(:,2)+0.5) <= dx/3 | abs(mesh_0.vtx(:,1)-1.5) <= dx/3 | abs(mesh_0.vtx(:,2)-1.5) <= dx/3;

intern_mask = mesh_0.vtx(:,1) >= 0 & mesh_0.vtx(:,1) <= 1 & mesh_0.vtx(:,2) >= 0 & mesh_0.vtx(:,2) <= 1;

u_til = u_list1(bound_mask,:);

a_exact = afun(mesh_0.vtx);

a_0 = ones(size(mesh_0.vtx(:,1)));


% %% Comparison between the post-processing methods
% 
% a_list_1 = gradient_1_final(u_til,a_0,gamma_0bis,qbis,0,0,N,mesh_0,dt,Nt,pfun,0,0,"None");
% a_list_2 = gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh_0,dt,Nt,pfun,0,"Smooth",1);
% a_list_3 = gradient_1_final(u_til,a_0,gamma_0,q,0,0,N,mesh_0,dt,Nt,pfun,0,0,"Cutoff");
% a_list_4 = gradient_1_final(u_til,a_0,gamma_0,q,0,0,N,mesh_0,dt,Nt,pfun,0,0,"Putup");
% 
% error_1 = sqrt(sum((a_list_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_2 = sqrt(sum((a_list_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3 = sqrt(sum((a_list_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4 = sqrt(sum((a_list_4(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
% % Plot the results
% figure
% scatter(0:N,error_1, "filled", DisplayName="None")
% hold on
% scatter(0:N,error_2, "filled", DisplayName="Smooth")
% scatter(0:N,error_3, "filled", DisplayName="Cut-off")
% scatter(0:N,error_4, "filled", DisplayName="Put up")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison of the post-processing methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% Vh = fem(mesh_0, 'P1');
% 
% figure
% graph(Vh, a_list_1(:,9))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 8$ (None)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_2(:,9))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 8$ (Smoothing)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_3(:,9))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 8$ (Cut-off)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_4(:,9))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 8$ (Put-up)","Interpreter","latex","FontSize",18)


% %% Simulation for the whole boundary
% 
% gamma_boundary_none = 0.3;
% gamma_boundary_cutoff = 0.04;
% gamma_boundary_putup = 0.08;
% qboundary = 1;
% 
% a_list_1 = gradient_1_final(u_til,a_0,gamma_boundary_none,qboundary,0,0,N,mesh_0,dt,Nt,pfun,1,0,"None");
% a_list_2 = gradient_1_1(u_til,a_0,gamma_boundary_none,qboundary,N,mesh_0,dt,Nt,pfun,1,"Smooth",1);
% a_list_3 = gradient_1_final(u_til,a_0,gamma_boundary_cutoff,qboundary,0,0,N,mesh_0,dt,Nt,pfun,1,0,"Cutoff");
% a_list_4 = gradient_1_final(u_til,a_0,gamma_boundary_putup,qboundary,0,0,N,mesh_0,dt,Nt,pfun,1,0,"Putup");
% 
% error_1 = sqrt(sum((a_list_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_2 = sqrt(sum((a_list_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3 = sqrt(sum((a_list_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4 = sqrt(sum((a_list_4(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
% % Plot the results
% figure
% scatter(0:N,error_1, "filled", DisplayName="None")
% hold on
% scatter(0:N,error_2, "filled", DisplayName="Smooth")
% scatter(0:N,error_3, "filled", DisplayName="Cut-off")
% scatter(0:N,error_4, "filled", DisplayName="Put up")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison of the post-processing methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% Vh = fem(mesh_0, 'P1');
% 
% figure
% graph(Vh, a_list_1(:,21))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 20$ (None)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_2(:,17))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 16$ (Smoothing)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_3(:,13))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 12$ (Cut-off)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_4(:,4))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 3$ (Put-up)","Interpreter","latex","FontSize",18)


% %% Simulation for the conjugate gradient
% 
% a_list_1 = conjugate_gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh_0,dt,Nt,pfun,0,"None",1);
% a_list_2 = conjugate_gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh_0,dt,Nt,pfun,0,"Smooth",1);
% a_list_3 = conjugate_gradient_1_1(u_til,a_0,gamma_0,q,N,mesh_0,dt,Nt,pfun,0,"Cutoff",1);
% a_list_4 = conjugate_gradient_1_1(u_til,a_0,gamma_0,q,N,mesh_0,dt,Nt,pfun,0,"Putup",1);
% 
% error_1 = sqrt(sum((a_list_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_2 = sqrt(sum((a_list_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3 = sqrt(sum((a_list_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4 = sqrt(sum((a_list_4(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
% % Plot the results
% figure
% scatter(0:N,error_1, "filled", DisplayName="None")
% hold on
% scatter(0:N,error_2, "filled", DisplayName="Smooth")
% scatter(0:N,error_3, "filled", DisplayName="Cut-off")
% scatter(0:N,error_4, "filled", DisplayName="Put up")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison of the post-processing methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% Vh = fem(mesh_0, 'P1');
% 
% figure
% graph(Vh, a_list_1(:,8))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 7$ (None)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_2(:,8))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 7$ (Smoothing)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_3(:,8))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 7$ (Cut-off)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_4(:,8))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 7$ (Put-up)","Interpreter","latex","FontSize",18)


% %% Simulation with a tolerance and mean
% 
% tol_putup = 0.6;
% 
% a_list_3 = gradient_1_final(u_til,a_0,gamma_0,q,0,0,N,mesh_0,dt,Nt,pfun,0,1,"Cutoff");
% a_list_4 = gradient_1_final(u_til,a_0,gamma_0,q,tol_putup,0,N,mesh_0,dt,Nt,pfun,0,1,"Putup");
% 
% error_3 = sqrt(sum((a_list_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4 = sqrt(sum((a_list_4(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
% % Plot the results
% figure
% scatter(0:N,error_3, "filled", DisplayName="Cut-off")
% hold on
% scatter(0:N,error_4, "filled", DisplayName="Put up")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison of the post-processing methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% Vh = fem(mesh_0, 'P1');
% 
% figure
% graph(Vh, a_list_3(:,13))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 12$ (Cut-off)","Interpreter","latex","FontSize",18)
% 
% figure
% graph(Vh, a_list_4(:,13))
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 12$ (Put-up)","Interpreter","latex","FontSize",18)


% %% Simulation for the adaptive gradient method
% 
% mesh_ref = mesh_0;
% 
% [u_list2, mesh_0, t_list] = forward_solver_method2_1(2^(-3),dt,Nt,afun,pfun,delta,theta);
% 
% bound_mask_bis = abs(mesh_0.vtx(:,2)-1.5) <= dx/3;
% 
% mesh_bound = mesh_0.vtx(bound_mask_bis,:);
% 
% u_til = zeros([nnz(bound_mask_bis),size(t_list,2)]);
% 
% for k = 1:size(u_til,1)
%     [~, ind] = min((mesh_ref.vtx(:,1)-mesh_bound(k,1)).^2.+(mesh_ref.vtx(:,2)-mesh_bound(k,2)).^2);
% 
%     u_til(k,:) = u_list1(ind,:);
% end
% 
% a_0 = ones(size(mesh_0.vtx(:,1)));
% 
% a_list_4 = adaptive_gradient_1_1(u_til,a_0,gamma_0,q,tol,N,mesh_0,dt,Nt,pfun,method,0,"Putup",1);
% 
% 
% % Plot the results
% Vh = fem(a_list_4{12,2}, 'P1');
% 
% figure
% graph(Vh, a_list_4{12,1})
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 11$ (Put-up)","Interpreter","latex","FontSize",18)
% 
% figure
% graph_mesh(a_list_4{12,1},a_list_4{12,2},1)
% colorbar()
% xlabel("$x$","Interpreter","latex","FontSize",18)
% ylabel("$y$","Interpreter","latex","FontSize",18)
% title("$a_{n,h}$ for $n = 11$ (Put up)","Interpreter","latex","FontSize",18)
