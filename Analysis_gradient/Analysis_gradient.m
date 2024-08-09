%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the parameters for the gradient method and the numerical
% techniques of reconstruction using boundary measurements.


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


%% Generation of the measurement for Model 1

% Solve the forward problem
[u_list1, mesh, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);

% Extraction of the values on a part of the boundary
step = mesh.stp;
dx = step(1);

% Only the top side
bound_mask = abs(mesh.vtx(:,2)-1.5) <= dx/3;

% The whole boundary for the related tests
% bound_mask = abs(mesh.vtx(:,1)+0.5) <= dx/3 | abs(mesh.vtx(:,2)+0.5) <= dx/3 | abs(mesh.vtx(:,1)-1.5) <= dx/3 | abs(mesh.vtx(:,2)-1.5) <= dx/3;

u_til = u_list1(bound_mask,:);

% Extraction of the values in Omega_1
intern_mask = mesh.vtx(:,1) >= 0 & mesh.vtx(:,1) <= 1 & mesh.vtx(:,2) >= 0 & mesh.vtx(:,2) <= 1;

% Generation of the initial value and the exact coefficient
a_exact = afun(mesh.vtx);

a_0 = ones(size(a_exact));


% %% Analysis of Gamma_0
% 
% % % Generation of the results
% % % Generation of the list of the values of gamma_0
% % gamma_list = cat(2,0.01:0.01:0.1,0.2:0.1:1);
% % 
% % 
% % error_list = zeros(size(gamma_list));
% % error_list_cut = zeros(size(gamma_list));
% % error_list_put = zeros(size(gamma_list));
% % 
% % for k = 1:size(gamma_list,2)
% %     a_list = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,0,"None",1);
% %     a_list_cut = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% %     a_list_put = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,0,"Putup",1);
% % 
% %     % Calculation of the relative quadratic error
% %     error_quad = sqrt(sum((a_list(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_cut = sqrt(sum((a_list_cut(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_put = sqrt(sum((a_list_put(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% %     % Selection of the minimal error over the iterations (does not include
% %     % the initial value)
% %     error_list(k) = min(error_quad(2:end));
% %     error_list_cut(k) = min(error_quad_cut(2:end));
% %     error_list_put(k) = min(error_quad_put(2:end));
% % end
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_gamma.mat 'gamma_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_gamma.mat 'gamma_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Plot the results
% figure
% scatter(gamma_list,error_list,"filled",DisplayName="None")
% hold on
% scatter(gamma_list,error_list_cut,"filled",DisplayName="Cut-off")
% scatter(gamma_list,error_list_put,"filled",DisplayName="Put-up")
% legend('Location','southeast')
% xlabel("$\gamma_0$","Interpreter","latex","FontSize",18)
% ylabel("$\min \frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of $\gamma_0$ for various post-processings","Interpreter","latex","FontSize",18)
% hold off


% %% Analysis of Gamma_0 (whole boundary)
% 
% % % Generation of the results
% % % Generation of the list of the values of gamma_0
% % gamma_list = cat(2,0.01:0.01:0.1,0.2:0.1:1);
% % 
% % 
% % error_list = zeros(size(gamma_list));
% % error_list_cut = zeros(size(gamma_list));
% % error_list_put = zeros(size(gamma_list));
% % 
% % for k = 1:size(gamma_list,2)
% %     a_list = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,1,"None",1);
% %     a_list_cut = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,1,"Cutoff",1);
% %     a_list_put = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,1,"Putup",1);
% % 
% %     % Calculation of the relative quadratic error
% %     error_quad = sqrt(sum((a_list(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_cut = sqrt(sum((a_list_cut(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_put = sqrt(sum((a_list_put(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% %     % Selection of the minimal error over the iterations (does not include
% %     % the initial value)
% %     error_list(k) = min(error_quad(2:end));
% %     error_list_cut(k) = min(error_quad_cut(2:end));
% %     error_list_put(k) = min(error_quad_put(2:end));
% % end
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_gamma_bound.mat 'gamma_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_gamma_bound.mat 'gamma_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Plot the results
% figure
% scatter(gamma_list,error_list,"filled",DisplayName="None")
% hold on
% scatter(gamma_list,error_list_cut,"filled",DisplayName="Cut-off")
% scatter(gamma_list,error_list_put,"filled",DisplayName="Put-up")
% legend('Location','southeast')
% xlabel("$\gamma_0$","Interpreter","latex","FontSize",18)
% ylabel("$\min \frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of $\gamma_0$ for various post-processings","Interpreter","latex","FontSize",18)
% hold off


% %% Analysis of q
% 
% % % Generation of the results
% % % Generation of the list of the values of q
% % q_list = 0.1:0.1:1;
% % 
% % 
% % error_list = zeros(size(q_list));
% % error_list_cut = zeros(size(q_list));
% % error_list_put = zeros(size(q_list));
% % 
% % for k = 1:size(gamma_list,2)
% %     a_list = gradient_1_1(u_til,a_0,gamma_0bis,q_list(k),N,mesh,dt,Nt,pfun,0,"None",1);
% %     a_list_cut = gradient_1_1(u_til,a_0,gamma_0,q_list(k),N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% %     a_list_put = gradient_1_1(u_til,a_0,gamma_0,q_list(k),N,mesh,dt,Nt,pfun,0,"Putup",1);
% % 
% %     % Calculation of the relative quadratic error
% %     error_quad = sqrt(sum((a_list(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_cut = sqrt(sum((a_list_cut(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_put = sqrt(sum((a_list_put(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% %     % Selection of the minimal error over the iterations (does not include
% %     % the initial value and the first iteration)
% %     error_list(k) = min(error_quad(3:end));
% %     error_list_cut(k) = min(error_quad_cut(3:end));
% %     error_list_put(k) = min(error_quad_put(3:end));
% % end
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_q.mat 'q_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_q.mat 'q_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Plot the results
% figure
% scatter(q_list,error_list,"filled",DisplayName="None")
% hold on
% scatter(q_list,error_list_cut,"filled",DisplayName="Cut-off")
% scatter(q_list,error_list_put,"filled",DisplayName="Put-up")
% legend('Location','southwest')
% xlabel("$q$","Interpreter","latex","FontSize",18)
% ylabel("$\min \frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of $q$ for various post-processings","Interpreter","latex","FontSize",18)
% hold off


% %% Analysis of q (whole boundary)
% 
% % % Generation of the results
% % % Generation of the list of the values of q
% % q_list = 0.1:0.1:1;
% % 
% % 
% % error_list = zeros(size(q_list));
% % error_list_cut = zeros(size(q_list));
% % error_list_put = zeros(size(q_list));
% % 
% % gamma_0 = 0.3;
% % gamma_0_cut = 0.04;
% % gamma_0_put = 0.08;
% % 
% % for k = 1:size(gamma_list,2)
% %     a_list = gradient_1_1(u_til,a_0,gamma_0,q_list(k),N,mesh,dt,Nt,pfun,1,"None",1);
% %     a_list_cut = gradient_1_1(u_til,a_0,gamma_0_cut,q_list(k),N,mesh,dt,Nt,pfun,1,"Cutoff",1);
% %     a_list_put = gradient_1_1(u_til,a_0,gamma_0_put,q_list(k),N,mesh,dt,Nt,pfun,1,"Putup",1);
% % 
% %     % Calculation of the relative quadratic error
% %     error_quad = sqrt(sum((a_list(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_cut = sqrt(sum((a_list_cut(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_put = sqrt(sum((a_list_put(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% %     % Selection of the minimal error over the iterations (does not include
% %     % the initial value and the first iteration)
% %     error_list(k) = min(error_quad(3:end));
% %     error_list_cut(k) = min(error_quad_cut(3:end));
% %     error_list_put(k) = min(error_quad_put(3:end));
% % end
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_q_bound.mat 'q_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_q_bound.mat 'q_list' 'error_list' 'error_list_cut' 'error_list_put'
% 
% % Plot the results
% figure
% scatter(q_list,error_list,"filled",DisplayName="None")
% hold on
% scatter(q_list,error_list_cut,"filled",DisplayName="Cut-off")
% scatter(q_list,error_list_put,"filled",DisplayName="Put-up")
% legend()
% xlabel("$q$","Interpreter","latex","FontSize",18)
% ylabel("$\min \frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of $q$ for various post-processings","Interpreter","latex","FontSize",18)
% hold off


% %% Comparison of the post-processing methods and the internal post-processing
% 
% % % Generation of the results
% % % Resolution for various post-processing methods
% % a_1_1 = gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh,dt,Nt,pfun,0,"None",0);
% % a_1_2 = gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh,dt,Nt,pfun,0,"None",1);
% % a_2_1 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Smooth",0);
% % a_2_2 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Smooth",1);
% % a_3_1 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",0);
% % a_3_2 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% % a_4_1 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",0);
% % a_4_2 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",1);
% % 
% % % Calculation of the error in Omega
% % error_1_1 = sqrt(sum((a_1_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_1_2 = sqrt(sum((a_1_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_2_1 = sqrt(sum((a_2_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_2_2 = sqrt(sum((a_2_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_3_1 = sqrt(sum((a_3_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_3_2 = sqrt(sum((a_3_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_4_1 = sqrt(sum((a_4_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_4_2 = sqrt(sum((a_4_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % 
% % % Calculation of the error in Omega_1
% % error_1_1bis = sqrt(sum((a_1_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_1_2bis = sqrt(sum((a_1_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_2_1bis = sqrt(sum((a_2_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_2_2bis = sqrt(sum((a_2_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_3_1bis = sqrt(sum((a_3_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_3_2bis = sqrt(sum((a_3_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_4_1bis = sqrt(sum((a_4_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_4_2bis = sqrt(sum((a_4_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_post_processing.mat error_1_1 error_1_2 error_2_1 error_2_2 error_3_1 error_3_2 ...
% %     error_4_1 error_4_2 error_1_1bis error_1_2bis error_2_1bis error_2_2bis error_3_1bis error_3_2bis error_4_1bis error_4_2bis
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_post_processing.mat error_1_1 error_1_2 error_2_1 error_2_2...
% error_3_1 error_3_2 error_4_1 error_4_2 error_1_1bis error_1_2bis error_2_1bis error_2_2bis error_3_1bis error_3_2bis error_4_1bis error_4_2bis
% 
% % Plot the results
% figure
% scatter(0:N,error_1_1, "filled", DisplayName="None, None")
% hold on
% scatter(0:N,error_1_2, LineWidth=0.9, DisplayName="None, Intern")
% scatter(0:N,error_2_1, "filled", DisplayName="Smooth, None")
% scatter(0:N,error_2_2, LineWidth=0.9, DisplayName="Smooth, Intern")
% scatter(0:N,error_3_1, "filled", DisplayName="Cut-off, None")
% scatter(0:N,error_3_2, LineWidth=0.9, DisplayName="Cut-off, Intern")
% scatter(0:N,error_4_1, "filled", DisplayName="Put up, None")
% scatter(0:N,error_4_2, LineWidth=0.9, DisplayName="Put up, Intern")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison of the post-processing methods and of the internal post-processing","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_1_1bis, "filled", DisplayName="None, None")
% hold on
% scatter(0:N,error_1_2bis, LineWidth=0.9, DisplayName="None, Intern")
% scatter(0:N,error_2_1bis, "filled", DisplayName="Smooth, None")
% scatter(0:N,error_2_2bis, LineWidth=0.9, DisplayName="Smooth, Intern")
% scatter(0:N,error_3_1bis, "filled", DisplayName="Cut-off, None")
% scatter(0:N,error_3_2bis, LineWidth=0.9, DisplayName="Cut-off, Intern")
% scatter(0:N,error_4_1bis, "filled", DisplayName="Put up, None")
% scatter(0:N,error_4_2bis, LineWidth=0.9, DisplayName="Put up, Intern")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison of the post-processing methods and of the internal post-processing","Interpreter","latex","FontSize",18)
% hold off


% %% Analysis of the integration time
% 
% % Generation of the list of the values of Nt
% Nt_list = 2000:500:10000;
% 
% % % Generation of the results
% % 
% % % Generation of the measurement
% % [u_list1, ~, ~] = forward_solver_method2_1(h,dt,Nt_list(end),afun,pfun,delta,theta);
% % 
% % error_1 = zeros(size(Nt_list));
% % error_3 = zeros(size(Nt_list));
% % error_4 = zeros(size(Nt_list));
% % 
% % error_1bis = zeros(size(Nt_list));
% % error_3bis = zeros(size(Nt_list));
% % error_4bis = zeros(size(Nt_list));
% % 
% % for k = 1:size(Nt_list,2)
% % 
% % u_til = u_list1(bound_mask,1:Nt_list(k)+1);
% % 
% % a_1 = gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh,dt,Nt_list(k),pfun,0,"None",1);
% % a_3 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt_list(k),pfun,0,"Cutoff",1);
% % a_4 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt_list(k),pfun,0,"Putup",1);
% % 
% % error_quad_1 = sqrt(sum((a_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_quad_3 = sqrt(sum((a_3 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_quad_4 = sqrt(sum((a_4 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % 
% % error_quad_1bis = sqrt(sum((a_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_quad_3bis = sqrt(sum((a_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_quad_4bis = sqrt(sum((a_4(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% % error_1(k) = min(error_quad_1);
% % error_3(k) = min(error_quad_3);
% % error_4(k) = min(error_quad_4);
% % 
% % error_1bis(k) = min(error_quad_1bis);
% % error_3bis(k) = min(error_quad_3bis);
% % error_4bis(k) = min(error_quad_4bis);
% % end
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_Nt.mat error_1 error_3 error_4 error_1bis error_3bis error_4bis
% 
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_Nt.mat error_1 error_3 error_4 error_1bis error_3bis error_4bis
% 
% % Plot the results
% figure
% scatter(Nt_list*dt,error_1, "filled", DisplayName="None")
% hold on
% scatter(Nt_list*dt,error_3, "filled", DisplayName="Cut-off")
% scatter(Nt_list*dt,error_4, "filled", DisplayName="Put up")
% legend('Location','southeast')
% xlabel("$T$","Interpreter","latex","FontSize",18)
% ylabel("$\min \frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of $T$ for various post-processings","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(Nt_list*dt,error_1bis, "filled", DisplayName="None")
% hold on
% scatter(Nt_list*dt,error_3bis, "filled", DisplayName="Cut-off")
% scatter(Nt_list*dt,error_4bis, "filled", DisplayName="Put up")
% legend('Location','southeast')
% xlabel("$T$","Interpreter","latex","FontSize",18)
% ylabel("$\min \frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of $T$ for various post-processings","Interpreter","latex","FontSize",18)
% hold off


% %% Analysis of the convergence in h
% 
% N = 10;
% 
% % % Generation of the results
% % % Generation of the measurement on a fine mesh
% % [u_list1, mesh_0, ~] = forward_solver_method2_1(2^(-6),dt,Nt,afun,pfun,delta,theta);
% % 
% % step = mesh_0.stp;
% % dx_0 = step(1);
% % 
% % bound_mask_0 = abs(mesh_0.vtx(:,2)-1.5) <= dx_0/3;
% % 
% % mesh_bound_0 = mesh_0.vtx(bound_mask_0,:);
% % 
% % u_til_1 = u_list1(bound_mask_0,:);
% % 
% % 
% % error_1 = zeros([5,N+1]);
% % error_3 = zeros([5,N+1]);
% % error_4 = zeros([5,N+1]);
% % error_1bis = zeros([5,N+1]);
% % error_3bis = zeros([5,N+1]);
% % error_4bis = zeros([5,N+1]);
% % 
% % for k = 2:6
% % 
% % % Creation of the mesh
% % NN = ceil(2*sqrt(2)/2^(-k));                                % Number of points, should be divisible by 4.
% % 
% % while mod(NN,4) ~= 0
% %     NN = NN+1;
% % end
% % 
% % mesh = mshSquare2(NN, [-0.5 1.5 -0.5 1.5]);
% % 
% % step = mesh.stp;
% % dx = step(1);
% % 
% % bound_mask = abs(mesh.vtx(:,2)-1.5) <= dx/3;
% % 
% % intern_mask = mesh.vtx(:,1) >= 0 & mesh.vtx(:,1) <= 1 & mesh.vtx(:,2) >= 0 & mesh.vtx(:,2) <= 1;
% % 
% % % Extraction of the data on a coarser mesh
% % u_til = zeros([nnz(bound_mask), size(u_list1,2)]);
% % 
% % mesh_bound = mesh.vtx(bound_mask,:);
% % 
% % for j = 1:size(u_til,1)
% % 
% %     [~,ind] = min((mesh_bound_0(:,1) - mesh_bound(j,1)).^2 + (mesh_bound_0(:,2) - mesh_bound(j,2)).^2);
% % 
% %     u_til(j,:) = u_til_1(ind,:);
% % 
% % end
% % 
% % a_exact = afun(mesh.vtx);
% % 
% % a_0 = ones(size(a_exact));
% % 
% % 
% % a_1 = gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh,dt,Nt,pfun,0,"None",1);
% % a_3 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% % a_4 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",1);
% % 
% % 
% % error_1(k-1,:) = sqrt(sum((a_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_3(k-1,:) = sqrt(sum((a_3 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % error_4(k-1,:) = sqrt(sum((a_4 - a_exact).^2))/sqrt(sum(a_exact.^2));
% % 
% % 
% % error_1bis(k-1,:) = sqrt(sum((a_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_3bis(k-1,:) = sqrt(sum((a_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_4bis(k-1,:) = sqrt(sum((a_4(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % end
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_hbis.mat error_1 error_3 error_4 error_1bis error_3bis error_4bis
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_hbis.mat error_1 error_3 error_4 error_1bis error_3bis error_4bis
% 
% % Plot the results
% figure
% scatter(0:N,error_1(1,:), "filled", DisplayName="h = 2^{-2}")
% hold on
% scatter(0:N,error_1(2,:), "filled", DisplayName="h = 2^{-3}")
% scatter(0:N,error_1(3,:), "filled", DisplayName="h = 2^{-4}")
% scatter(0:N,error_1(4,:), "filled", DisplayName="h = 2^{-5}")
% scatter(0:N,error_1(5,:), "filled", DisplayName="h = 2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the convergence in $h$ (None)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_1bis(1,:), "filled", DisplayName="h = 2^{-2}")
% hold on
% scatter(0:N,error_1bis(2,:), "filled", DisplayName="h = 2^{-3}")
% scatter(0:N,error_1bis(3,:), "filled", DisplayName="h = 2^{-4}")
% scatter(0:N,error_1bis(4,:), "filled", DisplayName="h = 2^{-5}")
% scatter(0:N,error_1bis(5,:), "filled", DisplayName="h = 2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the convergence in $h$ (None)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_3(1,:), "filled", DisplayName="h = 2^{-2}")
% hold on
% scatter(0:N,error_3(2,:), "filled", DisplayName="h = 2^{-3}")
% scatter(0:N,error_3(3,:), "filled", DisplayName="h = 2^{-4}")
% scatter(0:N,error_3(4,:), "filled", DisplayName="h = 2^{-5}")
% scatter(0:N,error_3(5,:), "filled", DisplayName="h = 2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the convergence in $h$ (Cut-off)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_3bis(1,:), "filled", DisplayName="h = 2^{-2}")
% hold on
% scatter(0:N,error_3bis(2,:), "filled", DisplayName="h = 2^{-3}")
% scatter(0:N,error_3bis(3,:), "filled", DisplayName="h = 2^{-4}")
% scatter(0:N,error_3bis(4,:), "filled", DisplayName="h = 2^{-5}")
% scatter(0:N,error_3bis(5,:), "filled", DisplayName="h = 2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the convergence in $h$ (Cut-off)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_4(1,:), "filled", DisplayName="h = 2^{-2}")
% hold on
% scatter(0:N,error_4(2,:), "filled", DisplayName="h = 2^{-3}")
% scatter(0:N,error_4(3,:), "filled", DisplayName="h = 2^{-4}")
% scatter(0:N,error_4(4,:), "filled", DisplayName="h = 2^{-5}")
% scatter(0:N,error_4(5,:), "filled", DisplayName="h = 2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the convergence in $h$ (Put-up)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_4bis(1,:), "filled", DisplayName="h = 2^{-2}")
% hold on
% scatter(0:N,error_4bis(2,:), "filled", DisplayName="h = 2^{-3}")
% scatter(0:N,error_4bis(3,:), "filled", DisplayName="h = 2^{-4}")
% scatter(0:N,error_4bis(4,:), "filled", DisplayName="h = 2^{-5}")
% scatter(0:N,error_4bis(5,:), "filled", DisplayName="h = 2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the convergence in $h$ (Put-up)","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% % Calculation of the order of convergence
% order_1bis = zeros([4,size(error_1,2)]);
% order_3bis = zeros([4,size(error_1,2)]);
% order_4bis = zeros([4,size(error_1,2)]);
% 
% for k = 1:4
%     order_1bis(k,:) = log2(error_1bis(k,:) ./ error_1bis(k+1,:));
%     order_3bis(k,:) = log2(error_3bis(k,:) ./ error_3bis(k+1,:));
%     order_4bis(k,:) = log2(error_4bis(k,:) ./ error_4bis(k+1,:));
% end
% 
% % Plot the results
% figure
% scatter(0:N,order_1bis(1,:), "filled", DisplayName="h = 2^{-2}/2^{-3}")
% hold on
% scatter(0:N,order_1bis(2,:), "filled", DisplayName="h = 2^{-3}/2^{-4}")
% scatter(0:N,order_1bis(3,:), "filled", DisplayName="h = 2^{-4}/2^{-5}")
% scatter(0:N,order_1bis(4,:), "filled", DisplayName="h = 2^{-5}/2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("Numerical order in $h$","Interpreter","latex","FontSize",18)
% title("Numerical order of convergence in $h$ (None)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,order_3bis(1,:), "filled", DisplayName="h = 2^{-2}/2^{-3}")
% hold on
% scatter(0:N,order_3bis(2,:), "filled", DisplayName="h = 2^{-3}/2^{-4}")
% scatter(0:N,order_3bis(3,:), "filled", DisplayName="h = 2^{-4}/2^{-5}")
% scatter(0:N,order_3bis(4,:), "filled", DisplayName="h = 2^{-5}/2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("Numerical order in $h$","Interpreter","latex","FontSize",18)
% title("Numerical order of convergence in $h$ (Cut-off)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,order_4bis(1,:), "filled", DisplayName="h = 2^{-2}/2^{-3}")
% hold on
% scatter(0:N,order_4bis(2,:), "filled", DisplayName="h = 2^{-3}/2^{-4}")
% scatter(0:N,order_4bis(3,:), "filled", DisplayName="h = 2^{-4}/2^{-5}")
% scatter(0:N,order_4bis(4,:), "filled", DisplayName="h = 2^{-5}/2^{-6}")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("Numerical order in $h$","Interpreter","latex","FontSize",18)
% title("Numerical order of convergence in $h$ (Put-up)","Interpreter","latex","FontSize",18)
% hold off


% %% Analysis of the tolerance parameter
% 
% % % Generation of the results
% % % Generation of the list of the values of tolerances
% % tol_list = 0:0.1:1;
% % 
% % error_list_cut_1 = zeros(size(tol_list));
% % error_list_put_1 = zeros(size(tol_list));
% % error_list_cut_2 = zeros(size(tol_list));
% % error_list_put_2 = zeros(size(tol_list));
% % 
% % for k = 1:size(tol_list,2)
% %     a_list_cut_1 = gradient_1_final(u_til,a_0,gamma_0,q,tol_list(k),0,N,mesh,dt,Nt,pfun,0,0,"Cutoff");
% %     a_list_put_1 = gradient_1_final(u_til,a_0,gamma_0,q,tol_list(k),0,N,mesh,dt,Nt,pfun,0,0,"Putup");
% %     a_list_cut_2 = gradient_1_final(u_til,a_0,gamma_0,q,tol_list(k),0,N,mesh,dt,Nt,pfun,0,1,"Cutoff");
% %     a_list_put_2 = gradient_1_final(u_til,a_0,gamma_0,q,tol_list(k),0,N,mesh,dt,Nt,pfun,0,1,"Putup");
% % 
% %     error_quad_cut_1 = sqrt(sum((a_list_cut_1(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_put_1 = sqrt(sum((a_list_put_1(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_cut_2 = sqrt(sum((a_list_cut_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% %     error_quad_put_2 = sqrt(sum((a_list_put_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% %     error_list_cut_1(k) = min(error_quad_cut_1(3:end));
% %     error_list_put_1(k) = min(error_quad_put_1(3:end));
% %     error_list_cut_2(k) = min(error_quad_cut_2(3:end));
% %     error_list_put_2(k) = min(error_quad_put_2(3:end));
% % end
% % 
% % % Save the data for future analysis
% % save \Data_Analysis_gradient_1_1_tol.mat 'tol_list' 'error_list_cut_1' 'error_list_cut_2' 'error_list_put_1' 'error_list_put_2'
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_tol.mat 'tol_list' 'error_list_cut_1' 'error_list_cut_2' 'error_list_put_1' 'error_list_put_2'
% 
% % Plot the results
% figure
% scatter(tol_list,error_list_cut_1,"filled",DisplayName="Cut-off, None")
% hold on
% scatter(tol_list,error_list_put_1,"filled",DisplayName="Put-up, None")
% scatter(tol_list,error_list_cut_2,"filled",DisplayName="Cut-off, Mean")
% scatter(tol_list,error_list_put_2,"filled",DisplayName="Put-up, Mean")
% legend('Location','northwest')
% xlabel("$tol$","Interpreter","latex","FontSize",18)
% ylabel("$\min \frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of $tol$ for various post-processings","Interpreter","latex","FontSize",18)
% hold off


% %% Efficiency of the mean post-processing
% 
% % % Generation of the results
% % tol_none = 0;
% % tol_cutoff = 0;
% % tol_putup_none = 0.2;
% % tol_putup_mean = 0.6;
% % 
% % a_1 = gradient_1_final(u_til,a_0,gamma_0bis,qbis,tol_none,0,N,mesh,dt,Nt,pfun,0,0,"None");
% % a_1_2 = gradient_1_final(u_til,a_0,gamma_0bis,qbis,tol_none,0,N,mesh,dt,Nt,pfun,0,1,"None");
% % a_3 = gradient_1_final(u_til,a_0,gamma_0,q,tol_cutoff,0,N,mesh,dt,Nt,pfun,0,0,"Cutoff");
% % a_3_2 = gradient_1_final(u_til,a_0,gamma_0,q,tol_cutoff,0,N,mesh,dt,Nt,pfun,0,1,"Cutoff");
% % a_4 = gradient_1_final(u_til,a_0,gamma_0,q,tol_putup_none,0,N,mesh,dt,Nt,pfun,0,0,"Putup");
% % a_4_2 = gradient_1_final(u_til,a_0,gamma_0,q,tol_putup_mean,0,N,mesh,dt,Nt,pfun,0,1,"Putup");
% % 
% % error_1 = sqrt(sum((a_1(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_1_2 = sqrt(sum((a_1_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_3 = sqrt(sum((a_3(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_3_2 = sqrt(sum((a_3_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_4 = sqrt(sum((a_4(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % error_4_2 = sqrt(sum((a_4_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% % 
% % % Save the data for future analysis
% % save C:\Users\loren\Dropbox\Stage\Algorithmes\Data_Analysis_gradient_1_1_mean.mat error_1 error_1_2 error_3 error_3_2 error_4 error_4_2
% 
% % Load the data
% load \Data_simulation_gradient\Data_Analysis_gradient_1_1_mean.mat error_1 error_1_2 error_3 error_3_2 error_4 error_4_2
% 
% % Plot the results
% figure
% scatter(0:N,error_1, "filled", DisplayName="No mean")
% hold on
% scatter(0:N,error_1_2, "filled", DisplayName="Mean")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the mean post-processing (None)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_3, "filled", DisplayName="No mean")
% hold on
% scatter(0:N,error_3_2, "filled", DisplayName="Mean")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the mean post-processing (Cut-off)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_4, "filled", DisplayName="No mean")
% hold on
% scatter(0:N,error_4_2, "filled", DisplayName="Mean")
% legend()
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the mean post-processing (Put-up)","Interpreter","latex","FontSize",18)
% hold off


% %% Comparison between the simple and the conjugate gradient
% 
% a_1 = gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh,dt,Nt,pfun,0,"None",1);
% a_1_2 = conjugate_gradient_1_1(u_til,a_0,gamma_0bis,qbis,N,mesh,dt,Nt,pfun,0,"None",1);
% a_3 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% a_3_2 = conjugate_gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% a_4 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",1);
% a_4_2 = conjugate_gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",1);
% 
% error_1 = sqrt(sum((a_1(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_1_2 = sqrt(sum((a_1_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3 = sqrt(sum((a_3(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3_2 = sqrt(sum((a_3_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4 = sqrt(sum((a_4(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4_2 = sqrt(sum((a_4_2(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
% % Plot the results
% figure
% scatter(0:N,error_1, "filled", DisplayName="Simple")
% hold on
% scatter(0:N,error_1_2, "filled", DisplayName="Conjugate")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the conjugate gradient (None)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_3, "filled", DisplayName="Simple")
% hold on
% scatter(0:N,error_3_2, "filled", DisplayName="Conjugate")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the conjugate gradient (Cut-off)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_4, "filled", DisplayName="Simple")
% hold on
% scatter(0:N,error_4_2, "filled", DisplayName="Conjugate")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Analysis of the conjugate gradient (Put-up)","Interpreter","latex","FontSize",18)
% hold off


% %% Comparison with the whole boundary
% 
% % Generation of the results
% % Generation of the measurements
% [u_list, mesh, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);
% 
% step = mesh.stp;
% dx = step(1);
% 
% bound_mask1 = abs(mesh.vtx(:,2)-1.5) <= dx/3;
% bound_mask2 = abs(mesh.vtx(:,1)+0.5) <= dx/3 | abs(mesh.vtx(:,2)+0.5) <= dx/3 | abs(mesh.vtx(:,1)-1.5) <= dx/3 | abs(mesh.vtx(:,2)-1.5) <= dx/3;
% 
% intern_mask = mesh.vtx(:,1) >= 0 & mesh.vtx(:,1) <= 1 & mesh.vtx(:,2) >= 0 & mesh.vtx(:,2) <= 1;
% 
% u_til1 = u_list(bound_mask1,:);
% u_til2 = u_list(bound_mask2,:);
% 
% a_exact = afun(mesh.vtx);
% 
% a_0 = ones(size(a_exact));
% 
% gamma_boundary_none = 0.3;
% gamma_boundary_cutoff = 0.04;
% gamma_boundary_putup = 0.08;
% qboundary = 1;
% 
% 
% a_1 = gradient_1_1(u_til1,a_0,gamma_0bis,qbis,N,mesh,dt,Nt,pfun,0,"None",1);
% a_1_2 = gradient_1_1(u_til2,a_0,gamma_boundary_none,qboundary,N,mesh,dt,Nt,pfun,1,"None",1);
% a_2 = gradient_1_1(u_til1,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% a_2_2 = gradient_1_1(u_til2,a_0,gamma_boundary_cutoff,qboundary,N,mesh,dt,Nt,pfun,1,"Cutoff",1);
% a_3 = gradient_1_1(u_til1,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",1);
% a_3_2 = gradient_1_1(u_til2,a_0,gamma_boundary_putup,qboundary,N,mesh,dt,Nt,pfun,1,"Putup",1);
% 
% error_1 = sqrt(sum((a_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_1_2 = sqrt(sum((a_1_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_2 = sqrt(sum((a_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_2_2 = sqrt(sum((a_2_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3 = sqrt(sum((a_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3_2 = sqrt(sum((a_3_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
% 
% % Plot the results
% figure
% scatter(0:N,error_1, "filled", DisplayName="Top")
% hold on
% scatter(0:N,error_1_2, "filled", DisplayName="Whole")
% legend('Location','southwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison with the whole boundary (None)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_2, "filled", DisplayName="Top")
% hold on
% scatter(0:N,error_2_2, "filled", DisplayName="Whole")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison with the whole boundary (Cut-off)","Interpreter","latex","FontSize",18)
% hold off
% 
% figure
% scatter(0:N,error_3, "filled", DisplayName="Top")
% hold on
% scatter(0:N,error_3_2, "filled", DisplayName="Whole")
% legend('Location','northwest')
% xlabel("Index of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{n,h} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
% title("Comparison with the whole boundary (Put-up)","Interpreter","latex","FontSize",18)
% hold off
