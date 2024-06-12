%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Analysis of the gradient method

%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions;


%% Definition of a and parameters

% Direct problem

h = 2^(-5);                                             % Mesh size

dt = 0.001;                                             % Time step (in seconds)

Nt = 2000;                                              % Number of time steps

om = 80;                                                % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

N = 15;

gamma_0 = 0.01;

q = 1;

tol = 0.95;

method = "linear";



%% Generation of the measurement

[u_list1, mesh, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,1/2,1/2);

step = mesh.stp;
dx = step(1);

% bound_mask = abs(mesh.vtx(:,1)+0.5) <= dx/3 | abs(mesh.vtx(:,2)+0.5) <= dx/3 | abs(mesh.vtx(:,1)-1.5) <= dx/3 | abs(mesh.vtx(:,2)-1.5) <= dx/3;
bound_mask = abs(mesh.vtx(:,2)-1.5) <= dx/3;

intern_mask = mesh.vtx(:,1) >= 0 & mesh.vtx(:,1) <= 1 & mesh.vtx(:,2) >= 0 & mesh.vtx(:,2) <= 1;

u_til = u_list1(bound_mask,:);

a_exact = afun(mesh.vtx);

a_0 = ones(size(a_exact));



% %% Gamma_0
% 
% gamma_list = cat(2,0.01:0.01:0.1,0.2:0.1:1);
% error_list = zeros(size(gamma_list));
% error_list_cut = zeros(size(gamma_list));
% error_list_put = zeros(size(gamma_list));
% 
% for k = 1:size(gamma_list,2)
%     a_list = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,0,"None",0);
%     a_list_cut = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,0,"Cutoff",0);
%     a_list_put = gradient_1_1(u_til,a_0,gamma_list(k),q,N,mesh,dt,Nt,pfun,0,"Putup",0);
% 
%     error_quad = sqrt(sum((a_list(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
%     error_quad_cut = sqrt(sum((a_list_cut(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
%     error_quad_put = sqrt(sum((a_list_put(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
%     error_list(k) = min(error_quad(2:end));
%     error_list_cut(k) = min(error_quad_cut(2:end));
%     error_list_put(k) = min(error_quad_put(2:end));
% 
%     disp(k);
% end
% 
% figure
% scatter(gamma_list,error_list,"filled",DisplayName="None")
% hold on
% scatter(gamma_list,error_list_cut,"filled",DisplayName="Cut-off")
% scatter(gamma_list,error_list_put,"filled",DisplayName="Put-up")
% legend()
% xlabel("$\gamma_0$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various $\gamma_0$","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% save C:\Users\loren\Dropbox\Stage\Algorithmes\Data_Analysis_gradient_1_1_gamma.mat 'gamma_list' 'error_list' 'error_list_cut' 'error_list_put'





% %% q
% 
% q_list = 0.1:0.1:1;
% error_list = zeros(size(q_list));
% error_list_cut = zeros(size(q_list));
% error_list_put = zeros(size(q_list));
% 
% for k = 1:size(q_list,2)
%     a_list = gradient_1_1(u_til,a_0,gamma_0,q_list(k),N,mesh,dt,Nt,pfun,0,"None",0);
%     a_list_cut = gradient_1_1(u_til,a_0,gamma_0,q_list(k),N,mesh,dt,Nt,pfun,0,"Cutoff",0);
%     a_list_put = gradient_1_1(u_til,a_0,gamma_0,q_list(k),N,mesh,dt,Nt,pfun,0,"Putup",0);
% 
%     error_quad = sqrt(sum((a_list(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
%     error_quad_cut = sqrt(sum((a_list_cut(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
%     error_quad_put = sqrt(sum((a_list_put(intern_mask,:)-a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
%     error_list(k) = min(error_quad(3:end));
%     error_list_cut(k) = min(error_quad_cut(3:end));
%     error_list_put(k) = min(error_quad_put(3:end));
% 
%     disp(k);
% end
% 
% figure
% scatter(q_list,error_list,"filled",DisplayName="None")
% hold on
% scatter(q_list,error_list_cut,"filled",DisplayName="Cut-off")
% scatter(q_list,error_list_put,"filled",DisplayName="Put-up")
% legend()
% xlabel("$q$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various $q$","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% save C:\Users\loren\Dropbox\Stage\Algorithmes\Data_Analysis_gradient_1_1_q.mat 'q_list' 'error_list' 'error_list_cut' 'error_list_put'





% %% Post_processing
% 
% N = 20;
% 
% a_1_1 = gradient_1_1(u_til,a_0,1,0.1,N,mesh,dt,Nt,pfun,0,"None",0);
% a_1_2 = gradient_1_1(u_til,a_0,1,0.1,N,mesh,dt,Nt,pfun,0,"None",1);
% a_2_1 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Smooth",0);
% a_2_2 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Smooth",1);
% a_3_1 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",0);
% a_3_2 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Cutoff",1);
% a_4_1 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",0);
% a_4_2 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt,pfun,0,"Putup",1);
% 
% 
% error_1_1 = sqrt(sum((a_1_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% error_1_2 = sqrt(sum((a_1_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% error_2_1 = sqrt(sum((a_2_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% error_2_2 = sqrt(sum((a_2_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% error_3_1 = sqrt(sum((a_3_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% error_3_2 = sqrt(sum((a_3_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% error_4_1 = sqrt(sum((a_4_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
% error_4_2 = sqrt(sum((a_4_2 - a_exact).^2))/sqrt(sum(a_exact.^2));
% 
% 
% error_1_1bis = sqrt(sum((a_1_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_1_2bis = sqrt(sum((a_1_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_2_1bis = sqrt(sum((a_2_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_2_2bis = sqrt(sum((a_2_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3_1bis = sqrt(sum((a_3_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_3_2bis = sqrt(sum((a_3_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4_1bis = sqrt(sum((a_4_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% error_4_2bis = sqrt(sum((a_4_2(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
% 
% 
% 
% 
% 
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
% xlabel("Rank of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various post-processing methods","Interpreter","latex","FontSize",18)
% hold off
% 
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
% xlabel("Rank of the iteration","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various post-processing methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% save C:\Users\loren\Dropbox\Stage\Algorithmes\Data_Analysis_gradient_1_1_post_processing.mat error_1_1 error_1_2 error_2_1 error_2_2...
% error_3_1 error_3_2 error_4_1 error_4_2 error_1_1bis error_1_2bis error_2_1bis error_2_2bis error_3_1bis error_3_2bis error_4_1bis error_4_2bis




%% Nt

Nt_list = 2000:500:10000;

error_1 = zeros(size(Nt_list));
error_3 = zeros(size(Nt_list));
error_4 = zeros(size(Nt_list));
error_1bis = zeros(size(Nt_list));
error_3bis = zeros(size(Nt_list));
error_4bis = zeros(size(Nt_list));

for k = 1:size(Nt_list,2)

[u_list1, ~, ~] = forward_solver_method2_1(h,dt,Nt_list(k),afun,pfun,1/2,1/2);

u_til = u_list1(bound_mask,:);

a_1 = gradient_1_1(u_til,a_0,1,0.1,N,mesh,dt,Nt_list(k),pfun,0,"None",1);
a_3 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt_list(k),pfun,0,"Cutoff",1);
a_4 = gradient_1_1(u_til,a_0,gamma_0,q,N,mesh,dt,Nt_list(k),pfun,0,"Putup",1);


error_quad_1 = sqrt(sum((a_1 - a_exact).^2))/sqrt(sum(a_exact.^2));
error_quad_3 = sqrt(sum((a_3 - a_exact).^2))/sqrt(sum(a_exact.^2));
error_quad_4 = sqrt(sum((a_4 - a_exact).^2))/sqrt(sum(a_exact.^2));


error_quad_1bis = sqrt(sum((a_1(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
error_quad_3bis = sqrt(sum((a_3(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));
error_quad_4bis = sqrt(sum((a_4(intern_mask,:) - a_exact(intern_mask)).^2))/sqrt(sum(a_exact(intern_mask).^2));

error_1(k) = min(error_quad_1);
error_3(k) = min(error_quad_3);
error_4(k) = min(error_quad_4);

error_1bis(k) = min(error_quad_1bis);
error_3bis(k) = min(error_quad_3bis);
error_4bis(k) = min(error_quad_4bis);

disp(k)
end





figure
scatter(Nt_list*dt,error_1, "filled", DisplayName="None")
hold on
scatter(Nt_list*dt,error_3, "filled", DisplayName="Cut-off")
scatter(Nt_list*dt,error_4, "filled", DisplayName="Put up")
legend()
xlabel("T","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
title("Relative quadratic error for various time integration","Interpreter","latex","FontSize",18)
hold off


figure
scatter(Nt_list*dt,error_1bis, "filled", DisplayName="None")
hold on
scatter(Nt_list*dt,error_3bis, "filled", DisplayName="Cut-off")
scatter(Nt_list*dt,error_4bis, "filled", DisplayName="Put up")
legend()
xlabel("T","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
title("Relative quadratic error for various time integration","Interpreter","latex","FontSize",18)
hold off



save C:\Users\loren\Dropbox\Stage\Algorithmes\Data_Analysis_gradient_1_1_Nt.mat error_1 error_3 error_4 error_1bis error_3bis error_4bis



