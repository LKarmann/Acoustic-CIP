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

dt = 0.001;                                             % Time step (in seconds)

Nt = 10000;                                             % Number of time steps

% t_list = 0:dt:Nt*dt;                                  % Time discretization

om = 80;                                                % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

delta = 1/2;                                            % Parameter of the numerical integration

theta = 1/2;                                            % Parameter of the numerical integration



%% Models for theta=1/2

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_4.mat t_list mesh U;
U1 = U;

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_8.mat U;
U2 = U;

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_12.mat U;
U3 = U;



% %% Comparison for various Nt
% 
% wh1 = discrete_laplace_transform(U1,t_list,10,"LeftRectangle");
% [Anh1, submesh] = inverse_solver(wh1, mesh, 10);
% A = afun(submesh.vtx);
% 
% Nt_list = 100:100:Nt;
% thres = 0.01;
% 
% 
% error_quad1 = zeros(size(Nt_list));
% error_quad2 = zeros(size(Nt_list));
% error_quad3 = zeros(size(Nt_list));
% 
% error_quad1_2 = zeros(size(Nt_list));
% error_quad2_2 = zeros(size(Nt_list));
% error_quad3_2 = zeros(size(Nt_list));
% 
% error_max1 = zeros(size(Nt_list));
% error_max2 = zeros(size(Nt_list));
% error_max3 = zeros(size(Nt_list));
% 
% error_max1_2 = zeros(size(Nt_list));
% error_max2_2 = zeros(size(Nt_list));
% error_max3_2 = zeros(size(Nt_list));
% 
% error_max_thres1 = zeros(size(Nt_list));
% error_max_thres2 = zeros(size(Nt_list));
% error_max_thres3 = zeros(size(Nt_list));
% 
% error_max_thres1_2 = zeros(size(Nt_list));
% error_max_thres2_2 = zeros(size(Nt_list));
% error_max_thres3_2 = zeros(size(Nt_list));
% 
% 
% for k = 1:size(Nt_list,2)
% 
%     wh1 = discrete_laplace_transform(U1(:,1:Nt_list(k)+1),0:dt:Nt_list(k)*dt,10,"LeftRectangle");
%     [Anh1, ~] = inverse_solver(wh1, mesh, 10);
%     Anh1_2 = boundary_post_processing(Anh1,submesh);
% 
%     wh2 = discrete_laplace_transform(U2(:,1:Nt_list(k)+1),0:dt:Nt_list(k)*dt,10,"LeftRectangle");
%     [Anh2, ~] = inverse_solver(wh2, mesh, 10);
%     Anh2_2 = boundary_post_processing(Anh2,submesh);
% 
%     wh3 = discrete_laplace_transform(U3(:,1:Nt_list(k)+1),0:dt:Nt_list(k)*dt,10,"LeftRectangle");
%     [Anh3, ~] = inverse_solver(wh3, mesh, 10);
%     Anh3_2 = boundary_post_processing(Anh3,submesh);
% 
% 
% 
%     error_quad1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
%     error_quad2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
%     error_quad3(k) = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% 
%     error_quad1_2(k) = sqrt(sum((Anh1_2-A).^2))/sqrt(sum(A.^2));
%     error_quad2_2(k) = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
%     error_quad3_2(k) = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
% 
% 
% 
%     error_max1(k) = max(abs(Anh1 - A));
%     error_max2(k) = max(abs(Anh2 - A));
%     error_max3(k) = max(abs(Anh3 - A));
% 
%     error_max1_2(k) = max(abs(Anh1_2 - A));
%     error_max2_2(k) = max(abs(Anh2_2 - A));
%     error_max3_2(k) = max(abs(Anh3_2 - A));
% 
% 
% 
%     B = unique(abs(Anh1 - A));
%     error_max_thres1(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh2 - A));
%     error_max_thres2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh3 - A));
%     error_max_thres3(k) = B(end-round(thres*size(B,1)));
% 
%     B = unique(abs(Anh1_2 - A));
%     error_max_thres1_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh2_2 - A));
%     error_max_thres2_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh3_2 - A));
%     error_max_thres3_2(k) = B(end-round(thres*size(B,1)));
% 
% 
% 
% end
% 
% 
% %% Comparison of errors
% 
% figure
% scatter(Nt_list*dt,error_quad1,'filled',DisplayName="Model 2.1")
% hold on
% scatter(Nt_list*dt,error_quad2,'filled',DisplayName="Model 2.2")
% scatter(Nt_list*dt,error_quad3,'filled',DisplayName="Model 2.3")
% axis([0,Nt*dt,0,1])
% legend()
% xlabel("T","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various integration time","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(Nt_list*dt,error_quad1_2,'filled',DisplayName="Model 2.1")
% hold on
% scatter(Nt_list*dt,error_quad2_2,'filled',DisplayName="Model 2.2")
% scatter(Nt_list*dt,error_quad3_2,'filled',DisplayName="Model 2.3")
% axis([0,Nt*dt,0,1])
% legend()
% xlabel("T","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error after post-processing for various integration time","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(Nt_list*dt,error_max1,'filled',DisplayName="Model 2.1")
% hold on
% scatter(Nt_list*dt,error_max2,'filled',DisplayName="Model 2.2")
% scatter(Nt_list*dt,error_max3,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("T","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error for various integration time","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(Nt_list*dt,error_max1_2,'filled',DisplayName="Model 2.1")
% hold on
% scatter(Nt_list*dt,error_max2_2,'filled',DisplayName="Model 2.2")
% scatter(Nt_list*dt,error_max3_2,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("T","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error after post-processing for various integration time","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(Nt_list*dt,error_max_thres1,'filled',DisplayName="Model 2.1")
% hold on
% scatter(Nt_list*dt,error_max_thres2,'filled',DisplayName="Model 2.2")
% scatter(Nt_list*dt,error_max_thres3,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("T","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold for various integration time","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(Nt_list*dt,error_max_thres1_2,'filled',DisplayName="Model 2.1")
% hold on
% scatter(Nt_list*dt,error_max_thres2_2,'filled',DisplayName="Model 2.2")
% scatter(Nt_list*dt,error_max_thres3_2,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("T","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold after post-processing for various integration time","Interpreter","latex","FontSize",18)
% hold off




% %% Comparison for various dt
% 
% wh1 = discrete_laplace_transform(U1,t_list,10,"LeftRectangle");
% [Anh1, submesh] = inverse_solver(wh1, mesh, 10);
% A = afun(submesh.vtx);
% 
% dt_list = dt:dt:100*dt;
% thres = 0.01;
% 
% 
% error_quad1 = zeros(size(dt_list));
% error_quad2 = zeros(size(dt_list));
% error_quad3 = zeros(size(dt_list));
% 
% error_quad1_2 = zeros(size(dt_list));
% error_quad2_2 = zeros(size(dt_list));
% error_quad3_2 = zeros(size(dt_list));
% 
% error_max1 = zeros(size(dt_list));
% error_max2 = zeros(size(dt_list));
% error_max3 = zeros(size(dt_list));
% 
% error_max1_2 = zeros(size(dt_list));
% error_max2_2 = zeros(size(dt_list));
% error_max3_2 = zeros(size(dt_list));
% 
% error_max_thres1 = zeros(size(dt_list));
% error_max_thres2 = zeros(size(dt_list));
% error_max_thres3 = zeros(size(dt_list));
% 
% error_max_thres1_2 = zeros(size(dt_list));
% error_max_thres2_2 = zeros(size(dt_list));
% error_max_thres3_2 = zeros(size(dt_list));
% 
% 
% for k = 1:size(dt_list,2)
% 
%     wh1 = discrete_laplace_transform(U1(:,1:k:end),t_list(1:k:end),10,"LeftRectangle");
%     [Anh1, ~] = inverse_solver(wh1, mesh, 10);
%     Anh1_2 = boundary_post_processing(Anh1,submesh);
% 
%     wh2 = discrete_laplace_transform(U2(:,1:k:end),t_list(1:k:end),10,"LeftRectangle");
%     [Anh2, ~] = inverse_solver(wh2, mesh, 10);
%     Anh2_2 = boundary_post_processing(Anh2,submesh);
% 
%     wh3 = discrete_laplace_transform(U3(:,1:k:end),t_list(1:k:end),10,"LeftRectangle");
%     [Anh3, ~] = inverse_solver(wh3, mesh, 10);
%     Anh3_2 = boundary_post_processing(Anh3,submesh);
% 
% 
% 
%     error_quad1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
%     error_quad2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
%     error_quad3(k) = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% 
%     error_quad1_2(k) = sqrt(sum((Anh1_2-A).^2))/sqrt(sum(A.^2));
%     error_quad2_2(k) = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
%     error_quad3_2(k) = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
% 
% 
% 
%     error_max1(k) = max(abs(Anh1 - A));
%     error_max2(k) = max(abs(Anh2 - A));
%     error_max3(k) = max(abs(Anh3 - A));
% 
%     error_max1_2(k) = max(abs(Anh1_2 - A));
%     error_max2_2(k) = max(abs(Anh2_2 - A));
%     error_max3_2(k) = max(abs(Anh3_2 - A));
% 
% 
% 
%     B = unique(abs(Anh1 - A));
%     error_max_thres1(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh2 - A));
%     error_max_thres2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh3 - A));
%     error_max_thres3(k) = B(end-round(thres*size(B,1)));
% 
%     B = unique(abs(Anh1_2 - A));
%     error_max_thres1_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh2_2 - A));
%     error_max_thres2_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh3_2 - A));
%     error_max_thres3_2(k) = B(end-round(thres*size(B,1)));
% 
% 
% 
% end
% 
% 
% %% Comparison of errors
% 
% figure
% scatter(dt_list,error_quad1,'filled',DisplayName="Model 2.1")
% hold on
% scatter(dt_list,error_quad2,'filled',DisplayName="Model 2.2")
% scatter(dt_list,error_quad3,'filled',DisplayName="Model 2.3")
% ylim([0,1])
% legend()
% xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various integration time step","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(dt_list,error_quad1_2,'filled',DisplayName="Model 2.1")
% hold on
% scatter(dt_list,error_quad2_2,'filled',DisplayName="Model 2.2")
% scatter(dt_list,error_quad3_2,'filled',DisplayName="Model 2.3")
% ylim([0,1])
% legend()
% xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error after post-processing for various integration time step","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(dt_list,error_max1,'filled',DisplayName="Model 2.1")
% hold on
% scatter(dt_list,error_max2,'filled',DisplayName="Model 2.2")
% scatter(dt_list,error_max3,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error for various integration time step","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(dt_list,error_max1_2,'filled',DisplayName="Model 2.1")
% hold on
% scatter(dt_list,error_max2_2,'filled',DisplayName="Model 2.2")
% scatter(dt_list,error_max3_2,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error after post-processing for various integration time step","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(dt_list,error_max_thres1,'filled',DisplayName="Model 2.1")
% hold on
% scatter(dt_list,error_max_thres2,'filled',DisplayName="Model 2.2")
% scatter(dt_list,error_max_thres3,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold for various integration time step","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(dt_list,error_max_thres1_2,'filled',DisplayName="Model 2.1")
% hold on
% scatter(dt_list,error_max_thres2_2,'filled',DisplayName="Model 2.2")
% scatter(dt_list,error_max_thres3_2,'filled',DisplayName="Model 2.3")
% legend()
% xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold after post-processing for various integration time step","Interpreter","latex","FontSize",18)
% hold off





% %% Comparison for various quadrature
% 
% wh1 = discrete_laplace_transform(U1,t_list,10,"LeftRectangle");
% [Anh1, submesh] = inverse_solver(wh1, mesh, 10);
% A = afun(submesh.vtx);
% 
% s = 0.5:0.5:30;
% thres = 0.01;
% 
% 
% error_quad1 = zeros(size(s));
% error_quad2 = zeros(size(s));
% error_quad3 = zeros(size(s));
% 
% error_quad1_2 = zeros(size(s));
% error_quad2_2 = zeros(size(s));
% error_quad3_2 = zeros(size(s));
% 
% error_max1 = zeros(size(s));
% error_max2 = zeros(size(s));
% error_max3 = zeros(size(s));
% 
% error_max1_2 = zeros(size(s));
% error_max2_2 = zeros(size(s));
% error_max3_2 = zeros(size(s));
% 
% error_max_thres1 = zeros(size(s));
% error_max_thres2 = zeros(size(s));
% error_max_thres3 = zeros(size(s));
% 
% error_max_thres1_2 = zeros(size(s));
% error_max_thres2_2 = zeros(size(s));
% error_max_thres3_2 = zeros(size(s));
% 
% error_quadbis1 = zeros(size(s));
% error_quadbis2 = zeros(size(s));
% error_quadbis3 = zeros(size(s));
% 
% error_quadbis1_2 = zeros(size(s));
% error_quadbis2_2 = zeros(size(s));
% error_quadbis3_2 = zeros(size(s));
% 
% error_maxbis1 = zeros(size(s));
% error_maxbis2 = zeros(size(s));
% error_maxbis3 = zeros(size(s));
% 
% error_maxbis1_2 = zeros(size(s));
% error_maxbis2_2 = zeros(size(s));
% error_maxbis3_2 = zeros(size(s));
% 
% error_max_thresbis1 = zeros(size(s));
% error_max_thresbis2 = zeros(size(s));
% error_max_thresbis3 = zeros(size(s));
% 
% error_max_thresbis1_2 = zeros(size(s));
% error_max_thresbis2_2 = zeros(size(s));
% error_max_thresbis3_2 = zeros(size(s));
% 
% error_quadter1 = zeros(size(s));
% error_quadter2 = zeros(size(s));
% error_quadter3 = zeros(size(s));
% 
% error_quadter1_2 = zeros(size(s));
% error_quadter2_2 = zeros(size(s));
% error_quadter3_2 = zeros(size(s));
% 
% error_maxter1 = zeros(size(s));
% error_maxter2 = zeros(size(s));
% error_maxter3 = zeros(size(s));
% 
% error_maxter1_2 = zeros(size(s));
% error_maxter2_2 = zeros(size(s));
% error_maxter3_2 = zeros(size(s));
% 
% error_max_threster1 = zeros(size(s));
% error_max_threster2 = zeros(size(s));
% error_max_threster3 = zeros(size(s));
% 
% error_max_threster1_2 = zeros(size(s));
% error_max_threster2_2 = zeros(size(s));
% error_max_threster3_2 = zeros(size(s));
% 
% 
% 
% for k = 1:size(s,2)
% 
%     wh1 = discrete_laplace_transform(U1,t_list,s(k),"LeftRectangle");
%     [Anh1, ~] = inverse_solver(wh1, mesh, s(k));
%     Anh1_2 = boundary_post_processing(Anh1,submesh);
% 
%     wh2 = discrete_laplace_transform(U2,t_list,s(k),"LeftRectangle");
%     [Anh2, ~] = inverse_solver(wh2, mesh, s(k));
%     Anh2_2 = boundary_post_processing(Anh2,submesh);
% 
%     wh3 = discrete_laplace_transform(U3,t_list,s(k),"LeftRectangle");
%     [Anh3, ~] = inverse_solver(wh3, mesh, s(k));
%     Anh3_2 = boundary_post_processing(Anh3,submesh);
% 
%     whbis1 = discrete_laplace_transform(U1,t_list,s(k),"RightRectangle");
%     [Anhbis1, ~] = inverse_solver(whbis1, mesh, s(k));
%     Anhbis1_2 = boundary_post_processing(Anhbis1,submesh);
% 
%     whbis2 = discrete_laplace_transform(U2,t_list,s(k),"RightRectangle");
%     [Anhbis2, ~] = inverse_solver(whbis2, mesh, s(k));
%     Anhbis2_2 = boundary_post_processing(Anhbis2,submesh);
% 
%     whbis3 = discrete_laplace_transform(U3,t_list,s(k),"RightRectangle");
%     [Anhbis3, ~] = inverse_solver(whbis3, mesh, s(k));
%     Anhbis3_2 = boundary_post_processing(Anhbis3,submesh);
% 
%     whter1 = discrete_laplace_transform(U1,t_list,s(k),"Trapezoidal");
%     [Anhter1, ~] = inverse_solver(whter1, mesh, s(k));
%     Anhter1_2 = boundary_post_processing(Anhter1,submesh);
% 
%     whter2 = discrete_laplace_transform(U2,t_list,s(k),"Trapezoidal");
%     [Anhter2, ~] = inverse_solver(whter2, mesh, s(k));
%     Anhter2_2 = boundary_post_processing(Anhter2,submesh);
% 
%     whter3 = discrete_laplace_transform(U3,t_list,s(k),"Trapezoidal");
%     [Anhter3, ~] = inverse_solver(whter3, mesh, s(k));
%     Anhter3_2 = boundary_post_processing(Anhter3,submesh);
% 
% 
% 
%     error_quad1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
%     error_quad2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
%     error_quad3(k) = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
% 
%     error_quad1_2(k) = sqrt(sum((Anh1_2-A).^2))/sqrt(sum(A.^2));
%     error_quad2_2(k) = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
%     error_quad3_2(k) = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
% 
% 
% 
%     error_max1(k) = max(abs(Anh1 - A));
%     error_max2(k) = max(abs(Anh2 - A));
%     error_max3(k) = max(abs(Anh3 - A));
% 
%     error_max1_2(k) = max(abs(Anh1_2 - A));
%     error_max2_2(k) = max(abs(Anh2_2 - A));
%     error_max3_2(k) = max(abs(Anh3_2 - A));
% 
% 
% 
%     B = unique(abs(Anh1 - A));
%     error_max_thres1(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh2 - A));
%     error_max_thres2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh3 - A));
%     error_max_thres3(k) = B(end-round(thres*size(B,1)));
% 
%     B = unique(abs(Anh1_2 - A));
%     error_max_thres1_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh2_2 - A));
%     error_max_thres2_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anh3_2 - A));
%     error_max_thres3_2(k) = B(end-round(thres*size(B,1)));
% 
% 
% 
% 
% 
%     error_quadbis1(k) = sqrt(sum((Anhbis1-A).^2))/sqrt(sum(A.^2));
%     error_quadbis2(k) = sqrt(sum((Anhbis2-A).^2))/sqrt(sum(A.^2));
%     error_quadbis3(k) = sqrt(sum((Anhbis3-A).^2))/sqrt(sum(A.^2));
% 
%     error_quadbis1_2(k) = sqrt(sum((Anhbis1_2-A).^2))/sqrt(sum(A.^2));
%     error_quadbis2_2(k) = sqrt(sum((Anhbis2_2-A).^2))/sqrt(sum(A.^2));
%     error_quadbis3_2(k) = sqrt(sum((Anhbis3_2-A).^2))/sqrt(sum(A.^2));
% 
% 
% 
%     error_maxbis1(k) = max(abs(Anhbis1 - A));
%     error_maxbis2(k) = max(abs(Anhbis2 - A));
%     error_maxbis3(k) = max(abs(Anhbis3 - A));
% 
%     error_maxbis1_2(k) = max(abs(Anhbis1_2 - A));
%     error_maxbis2_2(k) = max(abs(Anhbis2_2 - A));
%     error_maxbis3_2(k) = max(abs(Anhbis3_2 - A));
% 
% 
% 
%     B = unique(abs(Anhbis1 - A));
%     error_max_thresbis1(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhbis2 - A));
%     error_max_thresbis2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhbis3 - A));
%     error_max_thresbis3(k) = B(end-round(thres*size(B,1)));
% 
%     B = unique(abs(Anhbis1_2 - A));
%     error_max_thresbis1_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhbis2_2 - A));
%     error_max_thresbis2_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhbis3_2 - A));
%     error_max_thresbis3_2(k) = B(end-round(thres*size(B,1)));
% 
% 
% 
% 
% 
%     error_quadter1(k) = sqrt(sum((Anhter1-A).^2))/sqrt(sum(A.^2));
%     error_quadter2(k) = sqrt(sum((Anhter2-A).^2))/sqrt(sum(A.^2));
%     error_quadter3(k) = sqrt(sum((Anhter3-A).^2))/sqrt(sum(A.^2));
% 
%     error_quadter1_2(k) = sqrt(sum((Anhter1_2-A).^2))/sqrt(sum(A.^2));
%     error_quadter2_2(k) = sqrt(sum((Anhter2_2-A).^2))/sqrt(sum(A.^2));
%     error_quadter3_2(k) = sqrt(sum((Anhter3_2-A).^2))/sqrt(sum(A.^2));
% 
% 
% 
%     error_maxter1(k) = max(abs(Anhter1 - A));
%     error_maxter2(k) = max(abs(Anhter2 - A));
%     error_maxter3(k) = max(abs(Anhter3 - A));
% 
%     error_maxter1_2(k) = max(abs(Anhter1_2 - A));
%     error_maxter2_2(k) = max(abs(Anhter2_2 - A));
%     error_maxter3_2(k) = max(abs(Anhter3_2 - A));
% 
% 
% 
%     B = unique(abs(Anhter1 - A));
%     error_max_threster1(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhter2 - A));
%     error_max_threster2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhter3 - A));
%     error_max_threster3(k) = B(end-round(thres*size(B,1)));
% 
%     B = unique(abs(Anhter1_2 - A));
%     error_max_threster1_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhter2_2 - A));
%     error_max_threster2_2(k) = B(end-round(thres*size(B,1)));
%     B = unique(abs(Anhter3_2 - A));
%     error_max_threster3_2(k) = B(end-round(thres*size(B,1)));
% 
% 
% 
% end
% 
% 
% %% Comparison of errors
% 
% figure
% scatter(s,error_quad1,'filled',DisplayName="Left")
% hold on
% scatter(s,error_quadbis1,'filled',DisplayName="Right")
% scatter(s,error_quadter1,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_quad1_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_quadbis1_2,'filled',DisplayName="Right")
% scatter(s,error_quadter1_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max1,'filled',DisplayName="Left")
% hold on
% scatter(s,error_maxbis1,'filled',DisplayName="Right")
% scatter(s,error_maxter1,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max1_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_maxbis1_2,'filled',DisplayName="Right")
% scatter(s,error_maxter1_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max_thres1,'filled',DisplayName="Left")
% hold on
% scatter(s,error_max_thresbis1,'filled',DisplayName="Right")
% scatter(s,error_max_threster1,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max_thres1_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_max_thresbis1_2,'filled',DisplayName="Right")
% scatter(s,error_max_threster1_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% 
% 
% 
% figure
% scatter(s,error_quad2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_quadbis2,'filled',DisplayName="Right")
% scatter(s,error_quadter2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_quad2_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_quadbis2_2,'filled',DisplayName="Right")
% scatter(s,error_quadter2_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_maxbis2,'filled',DisplayName="Right")
% scatter(s,error_maxter2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max2_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_maxbis2_2,'filled',DisplayName="Right")
% scatter(s,error_maxter2_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max_thres2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_max_thresbis2,'filled',DisplayName="Right")
% scatter(s,error_max_threster2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max_thres2_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_max_thresbis2_2,'filled',DisplayName="Right")
% scatter(s,error_max_threster2_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% 
% 
% 
% figure
% scatter(s,error_quad3,'filled',DisplayName="Left")
% hold on
% scatter(s,error_quadbis3,'filled',DisplayName="Right")
% scatter(s,error_quadter3,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_quad3_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_quadbis3_2,'filled',DisplayName="Right")
% scatter(s,error_quadter3_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max3,'filled',DisplayName="Left")
% hold on
% scatter(s,error_maxbis3,'filled',DisplayName="Right")
% scatter(s,error_maxter3,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max3_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_maxbis3_2,'filled',DisplayName="Right")
% scatter(s,error_maxter3_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max_thres3,'filled',DisplayName="Left")
% hold on
% scatter(s,error_max_thresbis3,'filled',DisplayName="Right")
% scatter(s,error_max_threster3,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold for various quadrature","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max_thres3_2,'filled',DisplayName="Left")
% hold on
% scatter(s,error_max_thresbis3_2,'filled',DisplayName="Right")
% scatter(s,error_max_threster3_2,'filled',DisplayName="Trapezoidal")
% legend()
% xlabel("$s$","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold after post-processing for various quadrature","Interpreter","latex","FontSize",18)
% hold off




