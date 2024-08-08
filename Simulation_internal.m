%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Simulations and figures of the results from Forward solvers in time and
% reconstruction of the acoustic coefficient.


%% Preamble

clear;
close all;
format long;

addpath \gypsilab-master\gypsilab-master\openMsh;
addpath \gypsilab-master\gypsilab-master\openDom;
addpath \gypsilab-master\gypsilab-master\openFem;
addpath \Functions;


%% Definition of the parameters

% Forward problem

h = 2^(-5);                                             % Mesh size

dt = 0.001;                                             % Time step (in seconds)

Nt = 10000;                                             % Number of time steps

t_list = 0:dt:Nt*dt;                                    % Time discretization

om = 80;                                                % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source

% Test 1: 1 inclusion
afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% Test 2: 2 inclusions
% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

delta = 1/2;                                            % Parameter of the Newmark method

theta = 1/2;                                            % Parameter of the Newmark method


% Inverse problem

s = 5;                                                 % Pseudo-frequency


%% Solving the forward problem

[U_list2, mesh2, t_list2] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);

[U_list3, mesh3, t_list3] = forward_solver_method2_2(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);

[U_list4, mesh4, t_list4] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);


%% Plot the forward solution

% Construction of the Finite Elements
Vh = fem(mesh2,'P1');

figure
graph(Vh, U_list2(:,501))
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$u(x,t)$ for Test 1 at $t = 0.500$ $s$",'Interpreter','latex','FontSize',18)


% %% Movie of the solution
% 
% Vh = fem(mesh2,'P1');
% 
% figure
% hold on
% pause(3)
% for k = 1:20:size(U_list2,2)
%     clf
%     graph(Vh,U_list2(:,k))
%     colorbar()
%     drawnow
% end
% hold off


% %% Automatically save the figures
% 
% Vh = fem(mesh2,'P1');
% 
% z_min = min(U_list2(:));
% z_max = max(U_list2(:));
% 
% for n = 1:20:size(U_list2,2)
% 
% figure
% graph(Vh, U_list2(:,n));
% xlabel("$x$",'Interpreter','latex','FontSize',18)
% ylabel("$y$",'Interpreter','latex','FontSize',18)
% title("$u(x,t)$ at $t =$"+sprintf('%1.3f',(n-1)*dt)+" $s$",'Interpreter','latex','FontSize',18)
% clim([z_min, z_max])
% colorbar()
% 
% saveas(gcf,"Method_2_1_u_"+n+".png")
% close all
% end


%% Solving the inverse problem

% Discrete Laplace transform of the data
wh = discrete_laplace_transform(U_list4,t_list2,s,"Trapezoidal");

% Solve the inverse problem
[Ah, submesh2] = inverse_solver(wh,mesh2,s);

% Apply a post-processing on the reconstructed acoustic coefficient
Ah = boundary_post_processing(Ah,submesh2);


% Plot the result
Vh = fem(submesh2,'P1');

figure
graph(Vh, Ah)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ for Test 1 at $s=5$",'Interpreter','latex','FontSize',18)


%% Exact value of a

a_exact = afun(submesh2.vtx);

figure
graph(Vh,a_exact)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a(x,y)$ for Test 1",'Interpreter','latex','FontSize',18)


%% Calculation of the relative quadratic error and of the maximal error

error_quad = sqrt(sum((Ah-a_exact).^2))/sqrt(sum((a_exact).^2));

error_max = max(abs(Ah-a_exact));
