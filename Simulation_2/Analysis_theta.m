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



%% Model 2_1

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_1.mat t_list mesh U;
U1 = U;

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_2.mat U;
U2 = U;

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_3.mat U;
U3 = U;

load C:\Users\loren\Dropbox\Stage\Algorithmes\Simulation_2\Calculated_U\Data_U_Code_4.mat U;
U4 = U;



%% Calculation for various s

wh1 = discrete_laplace_transform(U1,t_list,10,"LeftRectangle");
[Anh1, submesh] = inverse_solver(wh1, mesh, 10);
A = afun(submesh.vtx);

s = 0.5:0.5:30;
thres = 0.01;


error_quad1 = zeros(size(s));
error_quad2 = zeros(size(s));
error_quad3 = zeros(size(s));
error_quad4 = zeros(size(s));

error_quad1_2 = zeros(size(s));
error_quad2_2 = zeros(size(s));
error_quad3_2 = zeros(size(s));
error_quad4_2 = zeros(size(s));

error_max1 = zeros(size(s));
error_max2 = zeros(size(s));
error_max3 = zeros(size(s));
error_max4 = zeros(size(s));

error_max1_2 = zeros(size(s));
error_max2_2 = zeros(size(s));
error_max3_2 = zeros(size(s));
error_max4_2 = zeros(size(s));

error_max_thres1 = zeros(size(s));
error_max_thres2 = zeros(size(s));
error_max_thres3 = zeros(size(s));
error_max_thres4 = zeros(size(s));

error_max_thres1_2 = zeros(size(s));
error_max_thres2_2 = zeros(size(s));
error_max_thres3_2 = zeros(size(s));
error_max_thres4_2 = zeros(size(s));


for k = 1:size(s,2)

    wh1 = discrete_laplace_transform(U1,t_list,s(k),"LeftRectangle");
    [Anh1, ~] = inverse_solver(wh1, mesh, s(k));
    Anh1_2 = boundary_post_processing(Anh1,submesh);

    wh2 = discrete_laplace_transform(U2,t_list,s(k),"LeftRectangle");
    [Anh2, ~] = inverse_solver(wh2, mesh, s(k));
    Anh2_2 = boundary_post_processing(Anh2,submesh);

    wh3 = discrete_laplace_transform(U3,t_list,s(k),"LeftRectangle");
    [Anh3, ~] = inverse_solver(wh3, mesh, s(k));
    Anh3_2 = boundary_post_processing(Anh3,submesh);

    wh4 = discrete_laplace_transform(U4,t_list,s(k),"LeftRectangle");
    [Anh4, ~] = inverse_solver(wh4, mesh, s(k));
    Anh4_2 = boundary_post_processing(Anh4,submesh);



    error_quad1(k) = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));
    error_quad2(k) = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));
    error_quad3(k) = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));
    error_quad4(k) = sqrt(sum((Anh4-A).^2))/sqrt(sum(A.^2));

    error_quad1_2(k) = sqrt(sum((Anh1_2-A).^2))/sqrt(sum(A.^2));
    error_quad2_2(k) = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));
    error_quad3_2(k) = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));
    error_quad4_2(k) = sqrt(sum((Anh4_2-A).^2))/sqrt(sum(A.^2));



    error_max1(k) = max(abs(Anh1 - A));
    error_max2(k) = max(abs(Anh2 - A));
    error_max3(k) = max(abs(Anh3 - A));
    error_max4(k) = max(abs(Anh4 - A));

    error_max1_2(k) = max(abs(Anh1_2 - A));
    error_max2_2(k) = max(abs(Anh2_2 - A));
    error_max3_2(k) = max(abs(Anh3_2 - A));
    error_max4_2(k) = max(abs(Anh4_2 - A));



    B = unique(abs(Anh1 - A));
    error_max_thres1(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh2 - A));
    error_max_thres2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh3 - A));
    error_max_thres3(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh4 - A));
    error_max_thres4(k) = B(end-round(thres*size(B,1)));

    B = unique(abs(Anh1_2 - A));
    error_max_thres1_2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh2_2 - A));
    error_max_thres2_2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh3_2 - A));
    error_max_thres3_2(k) = B(end-round(thres*size(B,1)));
    B = unique(abs(Anh4_2 - A));
    error_max_thres4_2(k) = B(end-round(thres*size(B,1)));



end





%% Comparison of errors

figure
scatter(s,error_quad1,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_quad2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_quad3,'filled',DisplayName="\theta = 1/4")
scatter(s,error_quad4,'filled',DisplayName="\theta = 1/2")
title("Relative quadratic error for various \theta")
hold off



figure
scatter(s,error_quad1_2,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_quad2_2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_quad3_2,'filled',DisplayName="\theta = 1/4")
scatter(s,error_quad4_2,'filled',DisplayName="\theta = 1/2")
title("Relative quadratic error after post-processing for various \theta")
hold off




figure
scatter(s,error_max1,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_max2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_max3,'filled',DisplayName="\theta = 1/4")
scatter(s,error_max4,'filled',DisplayName="\theta = 1/2")
title("Maximal error for various \theta")
hold off



figure
scatter(s,error_max1_2,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_max2_2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_max3_2,'filled',DisplayName="\theta = 1/4")
scatter(s,error_max4_2,'filled',DisplayName="\theta = 1/2")
title("Maximal error after post-processing for various \theta")
hold off




figure
scatter(s,error_max_thres1,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_max_thres2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_max_thres3,'filled',DisplayName="\theta = 1/4")
scatter(s,error_max_thres4,'filled',DisplayName="\theta = 1/2")
title("Maximal error with a threshold for various \theta")
hold off



figure
scatter(s,error_max_thres1_2,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_max_thres2_2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_max_thres3_2,'filled',DisplayName="\theta = 1/4")
scatter(s,error_max_thres4_2,'filled',DisplayName="\theta = 1/2")
title("Maximal error with a threshold after post-processing for various \theta")
hold off