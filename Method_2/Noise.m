%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Solves the Acoustic Coefficient Inverse Problem with Algorithm 2 and
% Method 2.

%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Fonctions;


%% Definition of a and parameters

% Direct problem
h = 2^(-5);                                             % Mesh size

om = 80;                                                % Frequency of the boudary source

noise_lvl = 3/100;                                      % Noise level

alph = 2*rand-1;                                        % Noise parameter

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);



% Inverse problem
s = 10;                                                % Pseudo-frequency

pfun1 = @(X) (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om) * (abs(X(:,2) - 1.5) < h/3);



%% Generation of the measurement for s_star and the relative error

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_4.mat t_list mesh U;
U1 = U;

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_12.mat U;
U2 = U;


[wh, ~] = forward_solver_method1(h, s, afun, pfun1);
wh1 = discrete_laplace_transform(U1,t_list,s,"Trapezoidal");
wh2 = discrete_laplace_transform(U2,t_list,s,"Trapezoidal");


[Anh_star, submesh] = inverse_solver(wh, mesh, s);
[Anh_star1, ~] = inverse_solver(wh1, mesh, s);
[Anh_star2, ~] = inverse_solver(wh2, mesh, s);

Anh_star_2 = boundary_post_processing(Anh_star,submesh);
Anh_star1_2 = boundary_post_processing(Anh_star1,submesh);
Anh_star2_2 = boundary_post_processing(Anh_star2,submesh);

Vh2 = fem(submesh, 'P1');

A = afun(submesh.vtx);




% %% Noise on wh
% 
% % %% Inversion without process
% % 
% % 
% % wh_til = wh*(1 + noise_lvl*alph);
% % wh_til1 = wh1*(1 + noise_lvl*alph);
% % wh_til2 = wh2*(1 + noise_lvl*alph);
% % 
% % [Anh, ~] = inverse_solver(wh_til, mesh, s);
% % [Anh1, ~] = inverse_solver(wh_til1, mesh, s);
% % [Anh2, ~] = inverse_solver(wh_til2, mesh, s);
% % 
% % 
% % figure
% % graph(Vh2, Anh)
% % title('Reconstruction of $a_{approx}$','interpreter','latex');
% % axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % figure
% % graph(Vh2, Anh1)
% % title('Reconstruction of $a_{approx}$','interpreter','latex');
% % axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % figure
% % graph(Vh2, Anh2)
% % title('Reconstruction of $a_{approx}$','interpreter','latex');
% % axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% 
% 
% 
% % %% Analysis of the error
% % 
% % alph_list = -1:0.01:1;
% % 
% % error_list = zeros(size(alph_list));
% % error_list1 = zeros(size(alph_list));
% % error_list2 = zeros(size(alph_list));
% % 
% % rel_error_list = zeros(size(alph_list));
% % rel_error_list1 = zeros(size(alph_list));
% % rel_error_list2 = zeros(size(alph_list));
% % 
% % 
% % for k = 1:size(alph_list,2)
% % 
% %     wh_til = wh*(1+noise_lvl*alph_list(k));
% %     wh_til1 = wh1*(1+noise_lvl*alph_list(k));
% %     wh_til2 = wh2*(1+noise_lvl*alph_list(k));
% % 
% % 
% %     [Anh, ~] = inverse_solver(wh_til, mesh, s);
% %     [Anh1, ~] = inverse_solver(wh_til1, mesh, s);
% %     [Anh2, ~] = inverse_solver(wh_til2, mesh, s);
% % 
% %     error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
% %     error_list1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
% %     error_list2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% % 
% %     rel_error_list(k) = sqrt(sum((Anh-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% %     rel_error_list1(k) = sqrt(sum((Anh1-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
% %     rel_error_list2(k) = sqrt(sum((Anh2-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));
% % 
% % end
% % 
% % 
% % 
% % 
% % figure
% % plot(alph_list, error_list)
% % hold on
% % yline(error_list(alph_list==0),'r--');
% % hold off
% % title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% % xlabel("\alpha","FontSize",15);
% % ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% % 
% % 
% % figure
% % plot(alph_list, error_list1)
% % hold on
% % yline(error_list1(alph_list==0),'r--');
% % hold off
% % title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% % xlabel("\alpha","FontSize",15);
% % ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% % 
% % 
% % figure
% % plot(alph_list, error_list2)
% % hold on
% % yline(error_list2(alph_list==0),'r--');
% % hold off
% % title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% % xlabel("\alpha","FontSize",15);
% % ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% % 
% % 
% % 
% % figure
% % plot(alph_list, rel_error_list)
% % hold on
% % plot(alph_list, rel_error_list1)
% % plot(alph_list, rel_error_list2)
% % title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% % xlabel("\alpha","FontSize",15);
% % ylabel("$\frac{||a_{approx} - a_{smooth}||_2}{||a_{smooth}||_2}$","Interpreter","latex","FontSize",18);
% 
% 
% 
% % %% Matrix noise uniformly distributed, V1
% % 
% % Alph = 2*rand(size(wh))-1;
% % 
% % wh_til = wh .* (1+noise_lvl*Alph);
% % 
% % 
% % 
% % % Without process
% % [Anh1, ~] = inverse_solver(wh_til, mesh, s);
% % 
% % 
% % error1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error1 = sqrt(sum((Anh1-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% % 
% % 
% % % 2D smooth
% % 
% % % 2D reconstruction
% % 
% % [~,~,vnh] = mesh2surface(wh_til,mesh);
% % 
% % 
% % % Smoothing method
% % vnh_smooth = smoothdata2(vnh,'gaussian');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % v = reshape(vnh_smooth,[],1);
% % 
% % 
% % % Solver
% % [Anh2, ~] = inverse_solver(v, mesh, s);
% % Anh2_2 = boundary_post_processing(Anh2,submesh);
% % 
% % 
% % error2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% % error2_2 = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error2 = sqrt(sum((Anh2-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% % rel_error2_2 = sqrt(sum((Anh2_2-Anh_star_2).^2))/sqrt(sum(Anh_star_2.^2));
% % 
% % 
% % 
% % 
% % 
% % figure
% % graph(Vh2,Anh_star)
% % title('Without noise','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % figure
% % graph(Vh2,Anh1)
% % title('Without pre-processing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % figure
% % graph(Vh2,Anh2)
% % title('With a 2D-gaussian smooth','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % 
% % % Second smoothing
% % 
% % [~,~,z] = mesh2surface(Anh2,submesh);
% % 
% % z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % Anh3 = reshape(z, [], 1);
% % 
% % 
% % error3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error3 = sqrt(sum((Anh3-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% % 
% % 
% % figure
% % graph(Vh2,Anh3)
% % title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % % Post-processing
% % 
% % 
% % figure
% % graph(Vh2,Anh2_2)
% % title('With a 2D-gaussian smooth','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % 
% % [~,~,z] = mesh2surface(Anh2_2,submesh);
% % 
% % z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % Anh3_2 = reshape(z, [], 1);
% % 
% % 
% % error3_2 = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error3_2 = sqrt(sum((Anh3_2-Anh_star_2).^2))/sqrt(sum(Anh_star_2.^2));
% % 
% % 
% % figure
% % graph(Vh2,Anh3_2)
% % title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% 
% 
% 
% 
% % %% Matrix noise uniformly distributed, V2
% % 
% % Alph = 2*rand(size(wh1))-1;
% % 
% % wh_til = wh1 .* (1+noise_lvl*Alph);
% % 
% % 
% % 
% % % Without process
% % [Anh1, ~] = inverse_solver(wh_til, mesh, s);
% % 
% % 
% % error1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error1 = sqrt(sum((Anh1-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
% % 
% % 
% % % 2D smooth
% % 
% % % 2D reconstruction
% % 
% % [~,~,vnh] = mesh2surface(wh_til,mesh);
% % 
% % 
% % % Smoothing method
% % vnh_smooth = smoothdata2(vnh,'gaussian');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % v = reshape(vnh_smooth,[],1);
% % 
% % 
% % % Solver
% % [Anh2, ~] = inverse_solver(v, mesh, s);
% % Anh2_2 = boundary_post_processing(Anh2,submesh);
% % 
% % 
% % error2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% % error2_2 = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error2 = sqrt(sum((Anh2-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
% % rel_error2_2 = sqrt(sum((Anh2_2-Anh_star1_2).^2))/sqrt(sum(Anh_star1_2.^2));
% % 
% % 
% % 
% % 
% % 
% % figure
% % graph(Vh2,Anh_star)
% % title('Without noise','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % figure
% % graph(Vh2,Anh1)
% % title('Without pre-processing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % figure
% % graph(Vh2,Anh2)
% % title('With a 2D-gaussian smooth','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % 
% % % Second smoothing
% % 
% % [~,~,z] = mesh2surface(Anh2,submesh);
% % 
% % z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % Anh3 = reshape(z, [], 1);
% % 
% % 
% % error3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error3 = sqrt(sum((Anh3-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
% % 
% % 
% % figure
% % graph(Vh2,Anh3)
% % title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % % Post-processing
% % 
% % 
% % figure
% % graph(Vh2,Anh2_2)
% % title('With a 2D-gaussian smooth','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % 
% % [~,~,z] = mesh2surface(Anh2_2,submesh);
% % 
% % z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % Anh3_2 = reshape(z, [], 1);
% % 
% % 
% % error3_2 = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error3_2 = sqrt(sum((Anh3_2-Anh_star1_2).^2))/sqrt(sum(Anh_star1_2.^2));
% % 
% % 
% % figure
% % graph(Vh2,Anh3_2)
% % title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% 
% 
% 
% 
% % %% Matrix noise uniformly distributed, V3
% % 
% % Alph = 2*rand(size(wh2))-1;
% % 
% % wh_til = wh2 .* (1+noise_lvl*Alph);
% % 
% % 
% % 
% % % Without process
% % [Anh1, ~] = inverse_solver(wh_til, mesh, s);
% % 
% % 
% % error1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error1 = sqrt(sum((Anh1-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));
% % 
% % 
% % % 2D smooth
% % 
% % % 2D reconstruction
% % 
% % [~,~,vnh] = mesh2surface(wh_til,mesh);
% % 
% % 
% % % Smoothing method
% % vnh_smooth = smoothdata2(vnh,'gaussian');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % v = reshape(vnh_smooth,[],1);
% % 
% % 
% % % Solver
% % [Anh2, ~] = inverse_solver(v, mesh, s);
% % Anh2_2 = boundary_post_processing(Anh2,submesh);
% % 
% % 
% % error2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% % error2_2 = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error2 = sqrt(sum((Anh2-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));
% % rel_error2_2 = sqrt(sum((Anh2_2-Anh_star2_2).^2))/sqrt(sum(Anh_star2_2.^2));
% % 
% % 
% % 
% % 
% % 
% % figure
% % graph(Vh2,Anh_star)
% % title('Without noise','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % figure
% % graph(Vh2,Anh1)
% % title('Without pre-processing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % figure
% % graph(Vh2,Anh2)
% % title('With a 2D-gaussian smooth','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % 
% % % Second smoothing
% % 
% % [~,~,z] = mesh2surface(Anh2,submesh);
% % 
% % z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % Anh3 = reshape(z, [], 1);
% % 
% % 
% % error3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error3 = sqrt(sum((Anh3-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));
% % 
% % 
% % figure
% % graph(Vh2,Anh3)
% % title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % % Post-processing
% % 
% % 
% % figure
% % graph(Vh2,Anh2_2)
% % title('With a 2D-gaussian smooth','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');
% % 
% % 
% % 
% % [~,~,z] = mesh2surface(Anh2_2,submesh);
% % 
% % z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% % 
% % Anh3_2 = reshape(z, [], 1);
% % 
% % 
% % error3_2 = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
% % 
% % rel_error3_2 = sqrt(sum((Anh3_2-Anh_star2_2).^2))/sqrt(sum(Anh_star2_2.^2));
% % 
% % 
% % figure
% % graph(Vh2,Anh3_2)
% % title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% % xlabel('$x$','interpreter','latex');
% % ylabel('$y$','interpreter','latex');
% % zlabel('$z$','interpreter','latex');





%% Noise on U

% %% Inversion without process
% 
% 
% wh_til = wh*(1 + noise_lvl*alph);
% U_til1 = U1*(1 + noise_lvl*alph);
% U_til2 = U2*(1 + noise_lvl*alph);
% 
% wh_til1 = discrete_laplace_transform(U_til1,t_list,s,"Trapezoidal");
% wh_til2 = discrete_laplace_transform(U_til2,t_list,s,"Trapezoidal");
% 
% [Anh, ~] = inverse_solver(wh_til, mesh, s);
% [Anh1, ~] = inverse_solver(wh_til1, mesh, s);
% [Anh2, ~] = inverse_solver(wh_til2, mesh, s);
% 
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
% 
% figure
% graph(Vh2, Anh2)
% title('Reconstruction of $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');



% %% Analysis of the error
% 
% alph_list = -1:0.01:1;
% 
% error_list = zeros(size(alph_list));
% error_list1 = zeros(size(alph_list));
% error_list2 = zeros(size(alph_list));
% 
% rel_error_list = zeros(size(alph_list));
% rel_error_list1 = zeros(size(alph_list));
% rel_error_list2 = zeros(size(alph_list));
% 
% 
% for k = 1:size(alph_list,2)
% 
%     wh_til = wh*(1+noise_lvl*alph_list(k));
%     U_til1 = U1*(1+noise_lvl*alph_list(k));
%     U_til2 = U2*(1+noise_lvl*alph_list(k));
% 
%     wh_til1 = discrete_laplace_transform(U_til1,t_list,s,"Trapezoidal");
%     wh_til2 = discrete_laplace_transform(U_til2,t_list,s,"Trapezoidal");
% 
% 
%     [Anh, ~] = inverse_solver(wh_til, mesh, s);
%     [Anh1, ~] = inverse_solver(wh_til1, mesh, s);
%     [Anh2, ~] = inverse_solver(wh_til2, mesh, s);
% 
%     error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
%     error_list1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
%     error_list2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% 
%     rel_error_list(k) = sqrt(sum((Anh-Anh_star).^2))/sqrt(sum(Anh_star.^2));
%     rel_error_list1(k) = sqrt(sum((Anh1-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
%     rel_error_list2(k) = sqrt(sum((Anh2-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));
% 
% end
% 
% 
% 
% 
% figure
% plot(alph_list, error_list)
% hold on
% yline(error_list(alph_list==0),'r--');
% hold off
% title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% xlabel("\alpha","FontSize",15);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% 
% 
% figure
% plot(alph_list, error_list1)
% hold on
% yline(error_list1(alph_list==0),'r--');
% hold off
% title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% xlabel("\alpha","FontSize",15);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% 
% 
% figure
% plot(alph_list, error_list2)
% hold on
% yline(error_list2(alph_list==0),'r--');
% hold off
% title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% xlabel("\alpha","FontSize",15);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% 
% 
% 
% figure
% plot(alph_list, rel_error_list)
% hold on
% plot(alph_list, rel_error_list1)
% plot(alph_list, rel_error_list2)
% title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% xlabel("\alpha","FontSize",15);
% ylabel("$\frac{||a_{approx} - a_{smooth}||_2}{||a_{smooth}||_2}$","Interpreter","latex","FontSize",18);



% %% Matrix noise uniformly distributed, V2
% 
% Alph = 2*rand(size(U1))-1;
% 
% U_til = U1 .* (1+noise_lvl*Alph);
% 
% wh_til = discrete_laplace_transform(U_til,t_list,s,"Trapezoidal");
% 
% 
% 
% % Without process
% [Anh1, ~] = inverse_solver(wh_til, mesh, s);
% 
% 
% error1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
% 
% rel_error1 = sqrt(sum((Anh1-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
% 
% 
% % 2D smooth
% 
% % 2D reconstruction
% 
% [~,~,vnh] = mesh2surface(wh_til,mesh);
% 
% 
% % Smoothing method
% vnh_smooth = smoothdata2(vnh,'gaussian');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% 
% v = reshape(vnh_smooth,[],1);
% 
% 
% % Solver
% [Anh2, ~] = inverse_solver(v, mesh, s);
% Anh2_2 = boundary_post_processing(Anh2,submesh);
% 
% 
% error2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% error2_2 = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
% 
% rel_error2 = sqrt(sum((Anh2-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
% rel_error2_2 = sqrt(sum((Anh2_2-Anh_star1_2).^2))/sqrt(sum(Anh_star1_2.^2));
% 
% 
% 
% 
% 
% figure
% graph(Vh2,Anh_star)
% title('Without noise','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2,Anh1)
% title('Without pre-processing','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2,Anh2)
% title('With a 2D-gaussian smooth','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% 
% 
% % Second smoothing
% 
% [~,~,z] = mesh2surface(Anh2,submesh);
% 
% z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% 
% Anh3 = reshape(z, [], 1);
% 
% 
% error3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% 
% rel_error3 = sqrt(sum((Anh3-Anh_star1).^2))/sqrt(sum(Anh_star1.^2));
% 
% 
% figure
% graph(Vh2,Anh3)
% title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% 
% % Post-processing
% 
% 
% figure
% graph(Vh2,Anh2_2)
% title('With a 2D-gaussian smooth','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% 
% 
% [~,~,z] = mesh2surface(Anh2_2,submesh);
% 
% z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% 
% Anh3_2 = reshape(z, [], 1);
% 
% 
% error3_2 = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
% 
% rel_error3_2 = sqrt(sum((Anh3_2-Anh_star1_2).^2))/sqrt(sum(Anh_star1_2.^2));
% 
% 
% figure
% graph(Vh2,Anh3_2)
% title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');



%% Matrix noise uniformly distributed, V3

Alph = 2*rand(size(U2))-1;

U_til = U2 .* (1+noise_lvl*Alph);

wh_til = discrete_laplace_transform(U_til,t_list,s,"Trapezoidal");



% Without process
[Anh1, ~] = inverse_solver(wh_til, mesh, s);


error1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));

rel_error1 = sqrt(sum((Anh1-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));


% 2D smooth

% 2D reconstruction

[~,~,vnh] = mesh2surface(wh_til,mesh);


% Smoothing method
vnh_smooth = smoothdata2(vnh,'gaussian');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'

v = reshape(vnh_smooth,[],1);


% Solver
[Anh2, ~] = inverse_solver(v, mesh, s);
Anh2_2 = boundary_post_processing(Anh2,submesh);


error2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
error2_2 = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));

rel_error2 = sqrt(sum((Anh2-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));
rel_error2_2 = sqrt(sum((Anh2_2-Anh_star2_2).^2))/sqrt(sum(Anh_star2_2.^2));





figure
graph(Vh2,Anh_star)
title('Without noise','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');

figure
graph(Vh2,Anh1)
title('Without pre-processing','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');

figure
graph(Vh2,Anh2)
title('With a 2D-gaussian smooth','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');



% Second smoothing

[~,~,z] = mesh2surface(Anh2,submesh);

z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'

Anh3 = reshape(z, [], 1);


error3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));

rel_error3 = sqrt(sum((Anh3-Anh_star2).^2))/sqrt(sum(Anh_star2.^2));


figure
graph(Vh2,Anh3)
title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');


% Post-processing


figure
graph(Vh2,Anh2_2)
title('With a 2D-gaussian smooth','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');



[~,~,z] = mesh2surface(Anh2_2,submesh);

z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'

Anh3_2 = reshape(z, [], 1);


error3_2 = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));

rel_error3_2 = sqrt(sum((Anh3_2-Anh_star2_2).^2))/sqrt(sum(Anh_star2_2.^2));


figure
graph(Vh2,Anh3_2)
title('With a pre-smoothing and a post-quadratic smoothing','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');











