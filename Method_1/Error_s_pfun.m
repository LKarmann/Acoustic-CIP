%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Analysis of the error on s and pfun.


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

s_star = 10.2;                                          % True pseudo-frequency

noise_lvl = 10/100;

s_inf = 0.5;

s_sup = 25;

s_step = 0.1;

tol = 0.001;

om = 80;                                                % Frequency of the boudary source

c = (1-exp(-2*pi*s_star/om))/((1+(s_star/om)^2)*om);    % Boundary source

Order = "increase";

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);


%% Generation of the measurement for s_star and the relative error

[wh, mesh] = forward_solver_method1(h, s_star, afun, pfun);



[Anh, submesh] = inverse_solver(wh, mesh, s_star);

A = afun(submesh.vtx);

error = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));



% %% Calculation of the error of reconstruction for various s
% 
% s_list = s_inf:s_step:s_sup;
% 
% error_list = zeros(size(s_list));
% 
% 
% for k = 1:size(s_list,2)
% 
%     s = s_list(k);
% 
%     [Anh, ~] = inverse_solver(wh, mesh, s);
% 
% 
%     error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
% 
% end
% 
% 
% 
% figure
% plot(s_list, error_list);
% hold on
% scatter(s_star,error,"filled");
% text(s_star,error+500,"$s^*$","Interpreter","latex");
% title("Relative quadratic error on $a$","Interpreter","latex","FontSize",18);
% xlabel("$s$","Interpreter","latex","FontSize",18);
% ylabel("$\frac{||a_{approx} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
% hold off
% 
% 
% figure
% plot(s_list, error_list);
% hold on
% scatter(s_star,error,"filled");
% axis([s_inf,s_sup,0,1]);
% text(s_star,error+0.02,"$s^*$","Interpreter","latex");
% xline(s_star*(1+noise_lvl),'-r');
% xline(s_star*(1-noise_lvl),'-r');
% title("Relative quadratic error on $a$ (zoom-in)","Interpreter","latex","FontSize",18);
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





% %% Blind reconstruction with an error on pfun
% 
% [~, ~, s_rec, error_rec] = blind_inverse_solver(wh,mesh,h,s_inf,s_sup,s_step,tol,pfun,Order);
% 
% 
% 
% c_noise = c*(1+noise_lvl*(2*rand-1));
% 
% pfun_noise = @(X) c_noise * (abs(X(:,2) - 1.5) < h/2);
% 
% [Anh, submesh, s_rec_2, error_rec_2] = blind_inverse_solver(wh,mesh,h,s_inf,s_sup,s_step,tol,pfun_noise,Order);
% 
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





% %% Further analysis on noise on c
% 
% c_list = linspace(c*(1-noise_lvl), c*(1+noise_lvl), 15);
% 
% error_list = zeros(size(c_list));
% 
% 
% figure
% hold on
% for k = 1:size(c_list,2)
% 
%     pfun = @(X) c_list(k) * (abs(X(:,2) - 1.5) < h/2);
% 
%     [Anh, submesh, s_rec_2, error_rec_2] = blind_inverse_solver(wh,mesh,h,s_inf,s_sup,s_step,tol,pfun,Order);
% 
%     scatter(c_list(k),s_rec_2,'blue','filled')
% 
%     error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
% 
% 
% end
% 
% title("Approximation of $s$ for different values of $c$","Interpreter","latex");
% xline(c,'r',"Exact c");
% xlabel("$c$",'Interpreter','latex');
% ylabel("$s$",'Interpreter','latex');
% hold off
% 
% 
% figure
% plot(c_list, error_list)
% title("Relative quadratic error on $a$","Interpreter","latex");
% xline(c,'r',"Exact c");
% xlabel("$c$",'Interpreter','latex');
% ylabel("$$\frac{||a_{approx} - a||_2}{||a||_2}$$",'Interpreter','latex');