%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the efficacy for various theta in the Newmark numerical 
% integration.


%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions;


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


%% Generation of the measurements in time for various theta

[U1, mesh, t_list] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,0);
[U2, ~, ~] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,1/12);
[U3, ~, ~] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,1/4);
[U4, ~, ~] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,1/2);


%% Reconstruction for various pseudo-frequencies

wh1 = discrete_laplace_transform(U1,t_list,10,"Trapezoidal");
[Ah1, submesh] = inverse_solver(wh1, mesh, 10);
a_exact = afun(submesh.vtx);

s = 1:0.1:30;

error_quad1 = zeros(size(s));
error_quad2 = zeros(size(s));
error_quad3 = zeros(size(s));
error_quad4 = zeros(size(s));

error_max1 = zeros(size(s));
error_max2 = zeros(size(s));
error_max3 = zeros(size(s));
error_max4 = zeros(size(s));


for k = 1:size(s,2)

    wh1 = discrete_laplace_transform(U1,t_list,s(k),"Trapezoidal");
    [Ah1, ~] = inverse_solver(wh1, mesh, s(k));
    Ah1 = boundary_post_processing(Ah1,submesh);

    wh2 = discrete_laplace_transform(U2,t_list,s(k),"Trapezoidal");
    [Ah2, ~] = inverse_solver(wh2, mesh, s(k));
    Ah2 = boundary_post_processing(Ah2,submesh);

    wh3 = discrete_laplace_transform(U3,t_list,s(k),"Trapezoidal");
    [Ah3, ~] = inverse_solver(wh3, mesh, s(k));
    Ah3 = boundary_post_processing(Ah3,submesh);

    wh4 = discrete_laplace_transform(U4,t_list,s(k),"Trapezoidal");
    [Ah4, ~] = inverse_solver(wh4, mesh, s(k));
    Ah4 = boundary_post_processing(Ah4,submesh);


    error_quad1(k) = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad2(k) = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad3(k) = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad4(k) = sqrt(sum((Ah4-a_exact).^2))/sqrt(sum(a_exact.^2));

    error_max1(k) = max(abs(Ah1 - a_exact));
    error_max2(k) = max(abs(Ah2 - a_exact));
    error_max3(k) = max(abs(Ah3 - a_exact));
    error_max4(k) = max(abs(Ah4 - a_exact));
end


%% Comparison between the thetas

figure
scatter(s,error_quad1,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_quad2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_quad3,'filled',DisplayName="\theta = 1/4")
scatter(s,error_quad4,'filled',DisplayName="\theta = 1/2")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Analysis of $\theta_{Newmark}$","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,error_quad1,'filled',DisplayName="\theta = 0")
hold on
scatter(s,error_quad2,'filled',DisplayName="\theta = 1/12")
scatter(s,error_quad3,'filled',DisplayName="\theta = 1/4")
scatter(s,error_quad4,'filled',DisplayName="\theta = 1/2")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Analysis of $\theta_{Newmark}$","Interpreter","latex","FontSize",18)
hold off


l1 = zeros(size(s));
l2 = zeros(size(s));

for k = 1:size(s,2)
    [~, ind1] = min([error_quad1(k), error_quad2(k), error_quad3(k), error_quad4(k)]);
    [~, ind2] = min([error_max1(k), error_max2(k), error_max3(k), error_max4(k)]);

    l1(k) = ind1;
    l2(k) = ind2;
end

figure
scatter(s,l1,'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,l2,'red',LineWidth=0.8,DisplayName="Maximal")
ylim([0.5,4.5])
yticks([1, 2, 3, 4])
yticklabels(["0", "1/12", "1/4", "1/2"])
xlabel("$s$","Interpreter","latex","FontSize",18)
legend()
title("Analysis of $\theta_{Newmark}$","Interpreter","latex","FontSize",18)
hold off