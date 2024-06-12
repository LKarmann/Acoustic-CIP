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
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions;


%% Definition of a and parameters

% Direct problem
h = 2^(-5);                                             % Mesh size

dt = 0.001;                                             % Time step (in seconds)

Nt = 10000;                                             % Number of time steps

% t_list = 0:dt:Nt*dt;                                  % Time discretization

om = 80;                                                % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source
% pfun1 = @(X) (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om) * (abs(X(:,2) - 1.5) < h/3);

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



%% Comparison between the forward solvers

wh1 = discrete_laplace_transform(U1,t_list,10,"LeftRectangle");
[Anh1, submesh] = inverse_solver(wh1, mesh, 10);
A = afun(submesh.vtx);

s = 1:0.1:30;
thres = 0.01;


error_quad = zeros(size(s));
error_quad1 = zeros(size(s));
error_quad2 = zeros(size(s));
error_quad3 = zeros(size(s));

error_quad_2 = zeros(size(s));
error_quad1_2 = zeros(size(s));
error_quad2_2 = zeros(size(s));
error_quad3_2 = zeros(size(s));

error_max = zeros(size(s));
error_max1 = zeros(size(s));
error_max2 = zeros(size(s));
error_max3 = zeros(size(s));

error_max_2 = zeros(size(s));
error_max1_2 = zeros(size(s));
error_max2_2 = zeros(size(s));
error_max3_2 = zeros(size(s));

error_max_thres = zeros(size(s));
error_max_thres1 = zeros(size(s));
error_max_thres2 = zeros(size(s));
error_max_thres3 = zeros(size(s));

error_max_thres_2 = zeros(size(s));
error_max_thres1_2 = zeros(size(s));
error_max_thres2_2 = zeros(size(s));
error_max_thres3_2 = zeros(size(s));


for k = 1:size(s,2)

    pfun1 = @(X) (1-exp(-2*pi*s(k)/om))/((1+(s(k)/om)^2)*om) * (abs(X(:,2) - 1.5) < h/3);

    [wh, ~] = forward_solver_method1(h,s(k),afun,pfun1);
    [Anh, ~] = inverse_solver(wh,mesh,s(k));
    Anh_2 = boundary_post_processing(Anh,submesh);

    wh1 = discrete_laplace_transform(U1,t_list,s(k),"Trapezoidal");
    [Anh1, ~] = inverse_solver(wh1, mesh, s(k));
    Anh1_2 = boundary_post_processing(Anh1,submesh);

    wh2 = discrete_laplace_transform(U2,t_list,s(k),"Trapezoidal");
    [Anh2, ~] = inverse_solver(wh2, mesh, s(k));
    Anh2_2 = boundary_post_processing(Anh2,submesh);

    wh3 = discrete_laplace_transform(U3,t_list,s(k),"Trapezoidal");
    [Anh3, ~] = inverse_solver(wh3, mesh, s(k));
    Anh3_2 = boundary_post_processing(Anh3,submesh);


    error_quad(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
    error_quad1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
    error_quad2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
    error_quad3(k) = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));

    error_quad_2(k) = sqrt(sum((Anh_2-A).^2))/sqrt(sum(A.^2));
    error_quad1_2(k) = sqrt(sum((Anh1_2-A).^2))/sqrt(sum(A.^2));
    error_quad2_2(k) = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
    error_quad3_2(k) = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));



    error_max(k) = max(abs(Anh - A));
    error_max1(k) = max(abs(Anh1 - A));
    error_max2(k) = max(abs(Anh2 - A));
    error_max3(k) = max(abs(Anh3 - A));

    error_max_2(k) = max(abs(Anh_2 - A));
    error_max1_2(k) = max(abs(Anh1_2 - A));
    error_max2_2(k) = max(abs(Anh2_2 - A));
    error_max3_2(k) = max(abs(Anh3_2 - A));



    B = unique(abs(Anh - A));
    error_max_thres(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh1 - A));
    error_max_thres1(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh2 - A));
    error_max_thres2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh3 - A));
    error_max_thres3(k) = B(end-round(thres*size(B,1)));

    B = unique(abs(Anh_2 - A));
    error_max_thres_2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh1_2 - A));
    error_max_thres1_2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh2_2 - A));
    error_max_thres2_2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh3_2 - A));
    error_max_thres3_2(k) = B(end-round(thres*size(B,1)));



end


% %% Comparison of errors
% 
% figure
% scatter(s,error_quad,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_quad1,'filled',DisplayName="Neumann")
% scatter(s,error_quad2,'filled',DisplayName="Mixed")
% scatter(s,error_quad3,'filled',DisplayName="Absorbing")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for various methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_quad_2,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_quad1_2,'filled',DisplayName="Neumann")
% scatter(s,error_quad2_2,'filled',DisplayName="Mixed")
% scatter(s,error_quad3_2,'filled',DisplayName="Absorbing")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error after post-processing for various methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max1,'filled',DisplayName="Neumann")
% scatter(s,error_max2,'filled',DisplayName="Mixed")
% scatter(s,error_max3,'filled',DisplayName="Absorbing")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error for various methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max_2,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max1_2,'filled',DisplayName="Neumann")
% scatter(s,error_max2_2,'filled',DisplayName="Mixed")
% scatter(s,error_max3_2,'filled',DisplayName="Absorbing")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error after post-processing for various methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max_thres,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max_thres1,'filled',DisplayName="Neumann")
% scatter(s,error_max_thres2,'filled',DisplayName="Mixed")
% scatter(s,error_max_thres3,'filled',DisplayName="Absorbing")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold for various methods","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max_thres_2,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max_thres1_2,'filled',DisplayName="Neumann")
% scatter(s,error_max_thres2_2,'filled',DisplayName="Mixed")
% scatter(s,error_max_thres3_2,'filled',DisplayName="Absorbing")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold after post-processing for various methods","Interpreter","latex","FontSize",18)
% hold off





% %% Method 1 vs 2.1
% 
% figure
% scatter(s,error_quad,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_quad1,'filled',DisplayName="Neumann")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error for Methods 1 and 2.2","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_quad_2,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_quad1_2,'filled',DisplayName="Neumann")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
% title("Relative quadratic error after post-processing for Methods 1 and 2.2","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max1,'filled',DisplayName="Neumann")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error for Methods 1 and 2.2","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max_2,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max1_2,'filled',DisplayName="Neumann")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
% title("Maximal error after post-processing for Methods 1 and 2.2","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% 
% figure
% scatter(s,error_max_thres,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max_thres1,'filled',DisplayName="Neumann")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold for Methods 1 and 2.2","Interpreter","latex","FontSize",18)
% hold off
% 
% 
% 
% figure
% scatter(s,error_max_thres_2,'filled',DisplayName="Laplace")
% hold on
% scatter(s,error_max_thres1_2,'filled',DisplayName="Neumann")
% legend()
% xlabel("s","Interpreter","latex","FontSize",18)
% ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
% title("Maximal error with a threshold after post-processing for Methods 1 and 2.2","Interpreter","latex","FontSize",18)
% hold off





%% Methods 2

figure
scatter(s,error_quad1,'filled',DisplayName="Neumann")
hold on
scatter(s,error_quad2,'filled',DisplayName="Mixed")
scatter(s,error_quad3,'filled',DisplayName="Absorbing")
legend()
xlabel("s","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
title("Relative quadratic error for various methods","Interpreter","latex","FontSize",18)
hold off



figure
scatter(s,error_quad1_2,'filled',DisplayName="Neumann")
hold on
scatter(s,error_quad2_2,'filled',DisplayName="Mixed")
scatter(s,error_quad3_2,'filled',DisplayName="Absorbing")
legend()
xlabel("s","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_{approx} - a_{exact}||_2}{||a_{exact}||_2}$","Interpreter","latex","FontSize",18)
title("Relative quadratic error after post-processing for various methods","Interpreter","latex","FontSize",18)
hold off




figure
scatter(s,error_max1,'filled',DisplayName="Neumann")
hold on
scatter(s,error_max2,'filled',DisplayName="Mixed")
scatter(s,error_max3,'filled',DisplayName="Absorbing")
legend()
xlabel("s","Interpreter","latex","FontSize",18)
ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
title("Maximal error for various methods","Interpreter","latex","FontSize",18)
hold off



figure
scatter(s,error_max1_2,'filled',DisplayName="Neumann")
hold on
scatter(s,error_max2_2,'filled',DisplayName="Mixed")
scatter(s,error_max3_2,'filled',DisplayName="Absorbing")
legend()
xlabel("s","Interpreter","latex","FontSize",18)
ylabel("$||a_{approx} - a_{exact}||_{\infty}$","Interpreter","latex","FontSize",18)
title("Maximal error after post-processing for various methods","Interpreter","latex","FontSize",18)
hold off




figure
scatter(s,error_max_thres1,'filled',DisplayName="Neumann")
hold on
scatter(s,error_max_thres2,'filled',DisplayName="Mixed")
scatter(s,error_max_thres3,'filled',DisplayName="Absorbing")
legend()
xlabel("s","Interpreter","latex","FontSize",18)
ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
title("Maximal error with a threshold for various methods","Interpreter","latex","FontSize",18)
hold off



figure
scatter(s,error_max_thres1_2,'filled',DisplayName="Neumann")
hold on
scatter(s,error_max_thres2_2,'filled',DisplayName="Mixed")
scatter(s,error_max_thres3_2,'filled',DisplayName="Absorbing")
legend()
xlabel("s","Interpreter","latex","FontSize",18)
ylabel("$||a_{approx} - a_{exact}||_{thres}$","Interpreter","latex","FontSize",18)
title("Maximal error with a threshold after post-processing for various methods","Interpreter","latex","FontSize",18)
hold off




figure
scatter(s,1+(error_quad1 > error_quad3),'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,1+(error_max1 > error_max3),'red',LineWidth=0.8,DisplayName="Maximal")
scatter(s,1+(error_max1_2 > error_max3_2),'+','green',LineWidth=0.8,DisplayName="Maximal post-processed")
ylim([0.5,2.5])
yticks([1, 2])
yticklabels(["Neumann", "Absorbing"])
xlabel("s","Interpreter","latex","FontSize",18)
legend()
title("Efficiency comparison between the forward methods","Interpreter","latex","FontSize",18)
hold off









