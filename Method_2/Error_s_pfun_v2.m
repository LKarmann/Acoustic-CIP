%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Analysis of the error on s and pfun with Method 2.


%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Fonctions\;


%% Definition of parameters

h = 2^(-5);                                             % Mesh size

dt = 0.001;

Nt = 10000;

s_star = 10.2;                                          % True pseudo-frequency

noise_lvl = 10/100;

s_inf = 5;

s_sup = 20;

s_step = 0.1;

tol = 0.001;

om = 80;                                                % Frequency of the boudary source

Order = "increase";

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

pfun = @(t) sin(om*t) .* (om * t <= 2*pi);               % Boundary source
pfun1 = @(X) (1-exp(-2*pi*s_star/om))/((1+(s_star/om)^2)*om) * (abs(X(:,2) - 1.5) < h/3);




%% Generation of the measurement for s_star and the relative error

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_4.mat t_list mesh U;
U1 = U;

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_12.mat U;
U2 = U;


[wh, ~] = forward_solver_method1(h, s_star, afun, pfun1);
wh1 = discrete_laplace_transform(U1,t_list,s_star,"Trapezoidal");
wh2 = discrete_laplace_transform(U2,t_list,s_star,"Trapezoidal");




[Anh, submesh] = inverse_solver(wh, mesh, s_star);
Anh_2 = boundary_post_processing(Anh,submesh);

A = afun(submesh.vtx);

error = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
error_2 = sqrt(sum((Anh_2-A).^2))/sqrt(sum(A.^2));



% %% Calculation of the error of reconstruction for various s
% 
% s_list = s_inf:s_step:s_sup;
% 
% error_list = zeros(size(s_list));
% error_list1 = zeros(size(s_list));
% error_list2 = zeros(size(s_list));
% error_list_2 = zeros(size(s_list));
% error_list1_2 = zeros(size(s_list));
% error_list2_2 = zeros(size(s_list));
% 
% 
% for k = 1:size(s_list,2)
% 
%     s = s_list(k);
% 
%     [Anh, ~] = inverse_solver(wh, mesh, s);
%     [Anh1, ~] = inverse_solver(wh1, mesh, s);
%     [Anh2, ~] = inverse_solver(wh2, mesh, s);
%     Anh_2 = boundary_post_processing(Anh,submesh);
%     Anh1_2 = boundary_post_processing(Anh1,submesh);
%     Anh2_2 = boundary_post_processing(Anh2,submesh);
% 
% 
%     error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
%     error_list1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
%     error_list2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% 
%     error_list_2(k) = sqrt(sum((Anh_2-A).^2))/sqrt(sum(A.^2));
%     error_list1_2(k) = sqrt(sum((Anh1_2-A).^2))/sqrt(sum(A.^2));
%     error_list2_2(k) = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
% 
% end
% 
% 
% 
% %% Plot the results
% 
% figure
% plot(s_list, error_list,DisplayName="Laplace");
% hold on
% plot(s_list, error_list1,DisplayName="Neuman");
% plot(s_list, error_list2,DisplayName="Absorbing");
% scatter(s_star,error,"filled");
% text(s_star,error+500,"$s^*$","Interpreter","latex");
% title("Relative quadratic error for various $s$","Interpreter","latex","FontSize",18);
% xlabel("$s$","Interpreter","latex","FontSize",18);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% hold off
% 
% 
% figure
% plot(s_list, error_list,DisplayName="Laplace");
% hold on
% plot(s_list, error_list1,DisplayName="Neuman");
% plot(s_list, error_list2,DisplayName="Absorbing");
% scatter(s_star,error,"filled");
% axis([s_inf,s_sup,0,1]);
% text(s_star,error+0.02,"$s^*$","Interpreter","latex");
% xline(s_star*(1+noise_lvl),'-r');
% xline(s_star*(1-noise_lvl),'-r');
% title("Relative quadratic error for various $s$ (zoom-in)","Interpreter","latex","FontSize",18);
% xlabel("$s$","Interpreter","latex","FontSize",18);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% hold off
% 
% 
% 
% figure
% plot(s_list, error_list_2,DisplayName="Laplace");
% hold on
% plot(s_list, error_list1_2,DisplayName="Neuman");
% plot(s_list, error_list2_2,DisplayName="Absorbing");
% scatter(s_star,error_2,"filled");
% text(s_star,error+500,"$s^*$","Interpreter","latex");
% title("Relative quadratic error for various $s$","Interpreter","latex","FontSize",18);
% xlabel("$s$","Interpreter","latex","FontSize",18);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% hold off
% 
% 
% figure
% plot(s_list, error_list_2,DisplayName="Laplace");
% hold on
% plot(s_list, error_list1_2,DisplayName="Neuman");
% plot(s_list, error_list2_2,DisplayName="Absorbing");
% scatter(s_star,error_2,"filled");
% axis([s_inf,s_sup,0,1]);
% text(s_star,error+0.02,"$s^*$","Interpreter","latex");
% xline(s_star*(1+noise_lvl),'-r');
% xline(s_star*(1-noise_lvl),'-r');
% title("Relative quadratic error for various $s$ (zoom-in)","Interpreter","latex","FontSize",18);
% xlabel("$s$","Interpreter","latex","FontSize",18);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% hold off
% 
% 
% 
% 
% s_noise = s_star*(1+noise_lvl*(2*rand -1));
% 
% [Anh, ~] = inverse_solver(wh, mesh, s_noise);
% [Anh1, ~] = inverse_solver(wh1, mesh, s_noise);
% [Anh2, ~] = inverse_solver(wh2, mesh, s_noise);
% 
% Vh2 = fem(submesh, 'P1');
% 
% figure
% graph(Vh2, Anh)
% title('Reconstruction of $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2, Anh1)
% title('Reconstruction of $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2, Anh2)
% title('Reconstruction of $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');



% %% Blind inverse solver
% 
% s_noise = s_star*(1+noise_lvl*(2*rand -1));
% pfun = @(t) sin(om * t) .* (om * t <= 2*pi);
% c = discrete_laplace_transform(pfun(t_list),t_list,s_noise,"Trapezoidal");
% pfunbis = @(X) c * (abs(X(:,2) - 1.5) < h/3);
% 
% [Anh, submesh, s, error_k] = blind_inverse_solver(wh2,mesh,h,s_inf,s_sup,s_step,tol,pfunbis,"increase");
% 
% 
% Vh2 = fem(submesh, 'P1');
% 
% figure
% graph(Vh2, Anh)
% title('Reconstruction of $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');



% %% Blind inverse solver with an error on pfun
% 
% s_list = s_star*(1+noise_lvl*(-1:0.1:1));
% 
% s_list1 = zeros(size(s_list));
% s_list2 = zeros(size(s_list));
% 
% error_list1 = zeros(size(s_list));
% error_list2 = zeros(size(s_list));
% 
% for k = 1:size(s_list,2)
% 
%     c = discrete_laplace_transform(pfun(t_list),t_list,s_list(k),"Trapezoidal");
%     pfunbis = @(X) c * (abs(X(:,2) - 1.5) < h/3);
% 
%     [Anh1, ~, s1, ~] = blind_inverse_solver(wh1,mesh,h,s_inf,s_sup,s_step,tol,pfunbis,"increase");
%     [Anh2, ~, s2, ~] = blind_inverse_solver(wh2,mesh,h,s_inf,s_sup,s_step,tol,pfunbis,"increase");
% 
%     s_list1(k) = s1;
%     s_list2(k) = s2;
% 
%     error_list1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
%     error_list2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% 
% end
% 
% 
% 
% 
% figure
% plot(s_list,s_list1)
% hold on
% plot(s_list,s_list2)
% plot(s_list,s_star*ones(size(s_list)),'k--')
% legend(["Neumann", "Absorbing", "s^*"])
% xlabel("$s$ for the calculation of pfun",'interpreter','latex',FontSize=18)
% ylabel("$s_{found}$",'interpreter','latex',FontSize=18)
% title("Blind inverse solver for various values of $s$ for the calculation of pfun",'interpreter','latex',FontSize=18)
% hold off
% 
% 
% figure
% plot(s_list,error_list1)
% hold on
% plot(s_list,error_list2)
% legend(["Neumann", "Absorbing"])
% xlabel("$s$ for the calculation of pfun",'interpreter','latex',FontSize=18)
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$",'interpreter','latex',FontSize=18)
% title("Blind inverse solver for various values of $s$ for the calculation of pfun",'interpreter','latex',FontSize=18)
% hold off






















