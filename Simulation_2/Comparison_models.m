%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the models and comparison between the forward solvers.


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


%% Generation of the measurements in time

[U1, mesh, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);
[U2, ~, ~] = forward_solver_method2_2(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);
[U3, ~, ~] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);



%% Reconstruction for various pseudo-frequencies

wh1 = discrete_laplace_transform(U1,t_list,10,"Trapezoidal");
[Ah1, submesh] = inverse_solver(wh1, mesh, 10);
a_exact = afun(submesh.vtx);

s = 1:0.1:30;

error_quad = zeros(size(s));
error_quad1 = zeros(size(s));
error_quad2 = zeros(size(s));
error_quad3 = zeros(size(s));

error_max = zeros(size(s));
error_max1 = zeros(size(s));
error_max2 = zeros(size(s));
error_max3 = zeros(size(s));


for k = 1:size(s,2)
    pfunbis = @(X) (1-exp(-2*pi*s(k)/om))/((1+(s(k)/om)^2)*om) * (abs(X(:,2) - 1.5) < h/3);

    [wh, ~] = forward_solver_method1(h,s(k),afun,pfunbis);
    [Ah, ~] = inverse_solver(wh,mesh,s(k));
    Ah = boundary_post_processing(Ah,submesh);

    wh1 = discrete_laplace_transform(U1,t_list,s(k),"Trapezoidal");
    [Ah1, ~] = inverse_solver(wh1, mesh, s(k));
    Ah1 = boundary_post_processing(Ah1,submesh);

    wh2 = discrete_laplace_transform(U2,t_list,s(k),"Trapezoidal");
    [Ah2, ~] = inverse_solver(wh2, mesh, s(k));
    Ah2 = boundary_post_processing(Ah2,submesh);

    wh3 = discrete_laplace_transform(U3,t_list,s(k),"Trapezoidal");
    [Ah3, ~] = inverse_solver(wh3, mesh, s(k));
    Ah3 = boundary_post_processing(Ah3,submesh);


    error_quad(k) = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad1(k) = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad2(k) = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad3(k) = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));

    error_max(k) = max(abs(Ah - a_exact));
    error_max1(k) = max(abs(Ah1 - a_exact));
    error_max2(k) = max(abs(Ah2 - a_exact));
    error_max3(k) = max(abs(Ah3 - a_exact));
end


%% Comparison between Method 1 and Forward solver 1

figure
scatter(s,error_quad,'filled',DisplayName="Laplace")
hold on
scatter(s,error_quad1,'filled',DisplayName="Neumann")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Comparison between Method 1 and Method 2","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,error_max,'filled',DisplayName="Laplace")
hold on
scatter(s,error_max1,'filled',DisplayName="Neumann")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Comparison between Method 1 and Method 2","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,1+(error_quad > error_quad1),'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,1+(error_max > error_max1),'red',LineWidth=0.8,DisplayName="Maximal")
ylim([0.5,2.5])
yticks([1, 2])
yticklabels(["Laplace", "Neumann"])
xlabel("$s$","Interpreter","latex","FontSize",18)
legend()
title("Efficiency comparison between the methods","Interpreter","latex","FontSize",18)
hold off


%% Comparison between the forward solvers in time

figure
scatter(s,error_quad1,'filled',DisplayName="Neumann")
hold on
scatter(s,error_quad2,'filled',DisplayName="Mixed")
scatter(s,error_quad3,'filled',DisplayName="Absorbing")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Comparison of the models","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,error_max1,'filled',DisplayName="Neumann")
hold on
scatter(s,error_max2,'filled',DisplayName="Mixed")
scatter(s,error_max3,'filled',DisplayName="Absorbing")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Comparison of the models","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,1+(error_quad1 > error_quad3),'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,1+(error_max1 > error_max3),'red',LineWidth=0.8,DisplayName="Maximal")
ylim([0.5,2.5])
yticks([1, 2])
yticklabels(["Neumann", "Absorbing"])
xlabel("$s$","Interpreter","latex","FontSize",18)
legend()
title("Efficiency comparison between the models","Interpreter","latex","FontSize",18)
hold off


l1 = zeros(size(s));
l2 = zeros(size(s));

for k = 1:size(s,2)
    [~, ind1] = min([error_quad1(k), error_quad2(k), error_quad3(k)]);
    [~, ind2] = min([error_max1(k), error_max2(k), error_max3(k)]);

    l1(k) = ind1;
    l2(k) = ind2;
end

figure
scatter(s,l1,'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,l2,'red',LineWidth=0.8,DisplayName="Maximal")
ylim([0.5,3.5])
yticks([1, 2, 3])
yticklabels(["Neumann", "Mixed", "Absorbing"])
xlabel("$s$","Interpreter","latex","FontSize",18)
legend()
title("Efficiency comparison between the models","Interpreter","latex","FontSize",18)
hold off
