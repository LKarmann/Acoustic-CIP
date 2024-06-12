%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Analysis of the noise on wh.


%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions\;


%% Definition of parameters

h = 2^(-5);                                             % Mesh size

s = 10;                                                 % True pseudo-frequency

noise_lvl = 3/100;                                      % Noise level

alph = 2*rand-1;                                        % Noise parameter

om = 80;                                                % Frequency of the boudary source

c = (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om);              % Boundary source

% afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);


%% Generation of the measurement

[wh, mesh] = forward_solver_method1(h, s, afun, pfun);


[Anh_star, submesh] = inverse_solver(wh, mesh, s);

A = afun(submesh.vtx);



% %% Inversion without process
% 
% 
% wh_til = wh*(1 + noise_lvl*alph);
% 
% [Anh, submesh] = inverse_solver(wh_til, mesh, s);
% 
% Vh2 = fem(submesh, 'P1');
% 
% 
% figure
% graph(Vh2, Anh)
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
% 
% rel_error_list = zeros(size(alph_list));
% 
% 
% for k = 1:size(alph_list,2)
% 
%     wh_til = wh*(1+noise_lvl*alph_list(k));
% 
% 
%     [Anh, ~] = inverse_solver(wh_til, mesh, s);
% 
%     error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
% 
%     rel_error_list(k) = sqrt(sum((Anh-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% end
% 
% figure
% plot(alph_list, error_list)
% yline(error_list(alph_list==0),'r');
% title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% xlabel("\alpha","FontSize",15);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% 
% 
% 
% figure
% plot(alph_list, rel_error_list)
% title("Relative quadratic error for noised measurements with \delta ="+sprintf('%1.0f', 100*noise_lvl)+" %");
% xlabel("\alpha","FontSize",15);
% ylabel("$\frac{||a_{approx} - a_{smooth}||_2}{||a_{smooth}||_2}$","Interpreter","latex","FontSize",18);



% %% Matrix noise uniformly distributed
% 
% Alph = 2*rand(size(wh))-1;
% 
% wh_til = wh .* (1+noise_lvl*Alph);
% 
% 
% 
% % Without process
% [Anh1, ~] = inverse_solver(wh_til, mesh, s);
% 
% 
% error1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
% 
% rel_error1 = sqrt(sum((Anh1-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% % 1D smooth
% [Anh2, ~] = inverse_solver(smoothdata(wh_til), mesh, s);
% 
% 
% error2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% 
% rel_error2 = sqrt(sum((Anh2-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% % 2D smooth
% 
% % 2D reconstruction
% x = unique(mesh.vtx(:,1));
% y = unique(mesh.vtx(:,2));
% 
% vnh = zeros(size(x,1));
% 
% for k = 1:size(vnh,1)
%     vnh(k,:) = wh_til(mesh.vtx(:,2)==y(k));
% end
% 
% 
% % Smoothing method
% vnh_smooth = smoothdata2(vnh,'gaussian');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% 
% 
% 
% % Solver
% [Anh3, ~] = inverse_solver(vnh_smooth, mesh, s);
% 
% 
% error3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% 
% rel_error3 = sqrt(sum((Anh3-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% 
% 
% 
% 
% Vh2 = fem(submesh,'P1');
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
% title('With a 1D smooth','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2,Anh3)
% title('With a 2D-gaussian smooth','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% 
% 
% % Second smoothing
% 
% [x,y,z] = submeshsurface(Anh3,submesh);
% 
% z = smoothdata2(z,'loess');         % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'
% 
% Anh4 = reshape(z, [], 1);
% 
% 
% error4 = sqrt(sum((Anh4-A).^2))/sqrt(sum(A.^2));
% 
% rel_error4 = sqrt(sum((Anh4-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% figure
% graph(Vh2,Anh4)
% title('With a pre-smoothing and a post-median smoothing','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');



% %% Matrix noise uniformly distributed
% 
% Alph = randn(size(wh));
% 
% wh_til = wh .* (1+noise_lvl*Alph);
% 
% 
% 
% % Without process
% [Anh1, ~] = inverse_solver(wh_til, mesh, s);
% 
% 
% error1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
% 
% rel_error1 = sqrt(sum((Anh1-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% % 1D smooth
% [Anh2, ~] = inverse_solver(smoothdata(wh_til), mesh, s);
% 
% 
% error2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
% 
% rel_error2 = sqrt(sum((Anh2-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% % 2D smooth
% 
% % 2D reconstruction
% x = unique(mesh.vtx(:,1));
% y = unique(mesh.vtx(:,2));
% 
% vnh = zeros(size(x,1));
% 
% for k = 1:size(vnh,1)
%     vnh(k,:) = wh_til(mesh.vtx(:,2)==y(k));
% end
% 
% 
% % Smoothing method
% vnh_smooth = smoothdata2(vnh,'gaussian');
% 
% 
% 
% % Solver
% [Anh3, ~] = inverse_solver(vnh_smooth, mesh, s);
% 
% 
% error3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% 
% rel_error3 = sqrt(sum((Anh3-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% 
% 
% 
% 
% Vh2 = fem(submesh,'P1');
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
% title('With a 1D smooth','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2,Anh3)
% title('With a 2D-gaussian smooth','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% 
% 
% % Second smoothing
% 
% [x,y,z] = submeshsurface(Anh3,submesh);
% 
% z = smoothdata2(z,'loess');
% 
% Anh4 = reshape(z, [], 1);
% 
% 
% error4 = sqrt(sum((Anh4-A).^2))/sqrt(sum(A.^2));
% 
% rel_error4 = sqrt(sum((Anh4-Anh_star).^2))/sqrt(sum(Anh_star.^2));
% 
% 
% figure
% graph(Vh2,Anh4)
% title('With a pre-smoothing and a post-median smoothing','interpreter','latex');
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');

