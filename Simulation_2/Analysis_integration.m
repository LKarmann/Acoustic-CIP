%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the numerical integration, the quadrature, the time step and
% the duration time.


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

% Inverse problem

s = 5;                                                 % Pseudo-frequency


%% Generation of the measurements in time

[U1, mesh, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);
[U2, ~, ~] = forward_solver_method2_2(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);
[U3, ~, ~] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);



%% Analysis of the quadrature

wh1 = discrete_laplace_transform(U1,t_list,10,"Trapezoidal");
[Ah1, submesh] = inverse_solver(wh1, mesh, 10);
a_exact = afun(submesh.vtx);

s = 1:0.1:30;

error_quad1 = zeros(size(s));
error_quad2 = zeros(size(s));
error_quad3 = zeros(size(s));

error_max1 = zeros(size(s));
error_max2 = zeros(size(s));
error_max3 = zeros(size(s));

error_quadbis1 = zeros(size(s));
error_quadbis2 = zeros(size(s));
error_quadbis3 = zeros(size(s));

error_maxbis1 = zeros(size(s));
error_maxbis2 = zeros(size(s));
error_maxbis3 = zeros(size(s));

error_quadter1 = zeros(size(s));
error_quadter2 = zeros(size(s));
error_quadter3 = zeros(size(s));

error_maxter1 = zeros(size(s));
error_maxter2 = zeros(size(s));
error_maxter3 = zeros(size(s));


for k = 1:size(s,2)

    wh1 = discrete_laplace_transform(U1,t_list,s(k),"LeftRectangle");
    [Ah1, ~] = inverse_solver(wh1, mesh, s(k));
    Ah1 = boundary_post_processing(Ah1,submesh);

    wh2 = discrete_laplace_transform(U2,t_list,s(k),"LeftRectangle");
    [Ah2, ~] = inverse_solver(wh2, mesh, s(k));
    Ah2 = boundary_post_processing(Ah2,submesh);

    wh3 = discrete_laplace_transform(U3,t_list,s(k),"LeftRectangle");
    [Ah3, ~] = inverse_solver(wh3, mesh, s(k));
    Ah3 = boundary_post_processing(Ah3,submesh);

    whbis1 = discrete_laplace_transform(U1,t_list,s(k),"RightRectangle");
    [Ahbis1, ~] = inverse_solver(whbis1, mesh, s(k));
    Ahbis1 = boundary_post_processing(Ahbis1,submesh);

    whbis2 = discrete_laplace_transform(U2,t_list,s(k),"RightRectangle");
    [Ahbis2, ~] = inverse_solver(whbis2, mesh, s(k));
    Ahbis2 = boundary_post_processing(Ahbis2,submesh);

    whbis3 = discrete_laplace_transform(U3,t_list,s(k),"RightRectangle");
    [Ahbis3, ~] = inverse_solver(whbis3, mesh, s(k));
    Ahbis3 = boundary_post_processing(Ahbis3,submesh);

    whter1 = discrete_laplace_transform(U1,t_list,s(k),"Trapezoidal");
    [Ahter1, ~] = inverse_solver(whter1, mesh, s(k));
    Ahter1 = boundary_post_processing(Ahter1,submesh);

    whter2 = discrete_laplace_transform(U2,t_list,s(k),"Trapezoidal");
    [Ahter2, ~] = inverse_solver(whter2, mesh, s(k));
    Ahter2 = boundary_post_processing(Ahter2,submesh);

    whter3 = discrete_laplace_transform(U3,t_list,s(k),"Trapezoidal");
    [Ahter3, ~] = inverse_solver(whter3, mesh, s(k));
    Ahter3 = boundary_post_processing(Ahter3,submesh);


    error_quad1(k) = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad2(k) = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad3(k) = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));

    error_max1(k) = max(abs(Ah1 - a_exact));
    error_max2(k) = max(abs(Ah2 - a_exact));
    error_max3(k) = max(abs(Ah3 - a_exact));

    error_quadbis1(k) = sqrt(sum((Ahbis1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quadbis2(k) = sqrt(sum((Ahbis2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quadbis3(k) = sqrt(sum((Ahbis3-a_exact).^2))/sqrt(sum(a_exact.^2));

    error_maxbis1(k) = max(abs(Ahbis1 - a_exact));
    error_maxbis2(k) = max(abs(Ahbis2 - a_exact));
    error_maxbis3(k) = max(abs(Ahbis3 - a_exact));

    error_quadter1(k) = sqrt(sum((Ahter1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quadter2(k) = sqrt(sum((Ahter2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quadter3(k) = sqrt(sum((Ahter3-a_exact).^2))/sqrt(sum(a_exact.^2));

    error_maxter1(k) = max(abs(Ahter1 - a_exact));
    error_maxter2(k) = max(abs(Ahter2 - a_exact));
    error_maxter3(k) = max(abs(Ahter3 - a_exact));
end


%% Plot the results

figure
scatter(s,error_quad1,'filled',DisplayName="Left")
hold on
scatter(s,error_quadbis1,'filled',DisplayName="Right")
scatter(s,error_quadter1,'filled',DisplayName="Trapezoidal")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Comparison of the quadrature","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,error_max1,'filled',DisplayName="Left")
hold on
scatter(s,error_maxbis1,'filled',DisplayName="Right")
scatter(s,error_maxter1,'filled',DisplayName="Trapezoidal")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Comparison of the quadrature","Interpreter","latex","FontSize",18)
hold off

l1 = zeros(size(s));
l2 = zeros(size(s));

for k = 1:size(s,2)
    [~, ind1] = min([error_quad1(k), error_quadbis1(k), error_quadter1(k)]);
    [~, ind2] = min([error_max1(k), error_maxbis1(k), error_maxter1(k)]);
    
    l1(k) = ind1;
    l2(k) = ind2;
end

figure
scatter(s,l1,'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,l2,'red',LineWidth=0.8,DisplayName="Maximal")
ylim([0.5,3.5])
yticks([1, 2, 3])
yticklabels(["Left rectangle", "Right rectangle", "Trapezoidal"])
xlabel("$s$","Interpreter","latex","FontSize",18)
legend()
title("Efficiency comparison between the quadrature rules","Interpreter","latex","FontSize",18)
hold off


figure
scatter(s,error_quad2,'filled',DisplayName="Left")
hold on
scatter(s,error_quadbis2,'filled',DisplayName="Right")
scatter(s,error_quadter2,'filled',DisplayName="Trapezoidal")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Comparison of the quadrature","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,error_max2,'filled',DisplayName="Left")
hold on
scatter(s,error_maxbis2,'filled',DisplayName="Right")
scatter(s,error_maxter2,'filled',DisplayName="Trapezoidal")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Comparison of the quadrature","Interpreter","latex","FontSize",18)
hold off

l1 = zeros(size(s));
l2 = zeros(size(s));

for k = 1:size(s,2)
    [~, ind1] = min([error_quad2(k), error_quadbis2(k), error_quadter2(k)]);
    [~, ind2] = min([error_max2(k), error_maxbis2(k), error_maxter2(k)]);
    
    l1(k) = ind1;
    l2(k) = ind2;
end

figure
scatter(s,l1,'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,l2,'red',LineWidth=0.8,DisplayName="Maximal")
ylim([0.5,3.5])
yticks([1, 2, 3])
yticklabels(["Left rectangle", "Right rectangle", "Trapezoidal"])
xlabel("$s$","Interpreter","latex","FontSize",18)
legend()
title("Efficiency comparison between the quadrature rules","Interpreter","latex","FontSize",18)
hold off


figure
scatter(s,error_quad3,'filled',DisplayName="Left")
hold on
scatter(s,error_quadbis3,'filled',DisplayName="Right")
scatter(s,error_quadter3,'filled',DisplayName="Trapezoidal")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Comparison of the quadrature","Interpreter","latex","FontSize",18)
hold off

figure
scatter(s,error_max3,'filled',DisplayName="Left")
hold on
scatter(s,error_maxbis3,'filled',DisplayName="Right")
scatter(s,error_maxter3,'filled',DisplayName="Trapezoidal")
legend()
xlabel("$s$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Comparison of the quadrature","Interpreter","latex","FontSize",18)
hold off

l1 = zeros(size(s));
l2 = zeros(size(s));

for k = 1:size(s,2)
    [~, ind1] = min([error_quad3(k), error_quadbis3(k), error_quadter3(k)]);
    [~, ind2] = min([error_max3(k), error_maxbis3(k), error_maxter3(k)]);
    
    l1(k) = ind1;
    l2(k) = ind2;
end

figure
scatter(s,l1,'filled','blue',DisplayName="Quadratic")
hold on
scatter(s,l2,'red',LineWidth=0.8,DisplayName="Maximal")
ylim([0.5,3.5])
yticks([1, 2, 3])
yticklabels(["Left rectangle", "Right rectangle", "Trapezoidal"])
xlabel("$s$","Interpreter","latex","FontSize",18)
legend()
title("Efficiency comparison between the quadrature rules","Interpreter","latex","FontSize",18)
hold off



%% Analysis of the time duration T

wh1 = discrete_laplace_transform(U1,t_list,10,"Trapezoidal");
[Ah1, submesh] = inverse_solver(wh1, mesh, 10);
a_exact = afun(submesh.vtx);

Nt_list = 100:100:Nt;

error_quad1 = zeros(size(Nt_list));
error_quad2 = zeros(size(Nt_list));
error_quad3 = zeros(size(Nt_list));

error_max1 = zeros(size(Nt_list));
error_max2 = zeros(size(Nt_list));
error_max3 = zeros(size(Nt_list));


for k = 1:size(Nt_list,2)

    wh1 = discrete_laplace_transform(U1(:,1:Nt_list(k)+1),0:dt:Nt_list(k)*dt,s,"Trapezoidal");
    [Ah1, ~] = inverse_solver(wh1, mesh, s);
    Ah1 = boundary_post_processing(Ah1,submesh);

    wh2 = discrete_laplace_transform(U2(:,1:Nt_list(k)+1),0:dt:Nt_list(k)*dt,s,"Trapezoidal");
    [Ah2, ~] = inverse_solver(wh2, mesh, s);
    Ah2 = boundary_post_processing(Ah2,submesh);

    wh3 = discrete_laplace_transform(U3(:,1:Nt_list(k)+1),0:dt:Nt_list(k)*dt,s,"Trapezoidal");
    [Ah3, ~] = inverse_solver(wh3, mesh, s);
    Ah3 = boundary_post_processing(Ah3,submesh);


    error_quad1(k) = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad2(k) = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad3(k) = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));

    error_max1(k) = max(abs(Ah1 - a_exact));
    error_max2(k) = max(abs(Ah2 - a_exact));
    error_max3(k) = max(abs(Ah3 - a_exact));
end


%% Plot the results

figure
scatter(Nt_list*dt,min(error_quad1,2),'filled',DisplayName="Neumann")
hold on
scatter(Nt_list*dt,min(error_quad2,2),'filled',DisplayName="Mixed")
scatter(Nt_list*dt,min(error_quad3,2),'filled',DisplayName="Absorbing")
legend()
xlabel("$T$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Error analysis for various $T$ (cut-off at 2)","Interpreter","latex","FontSize",18)
hold off

figure
scatter(Nt_list*dt,min(error_max1,3),'filled',DisplayName="Neumann")
hold on
scatter(Nt_list*dt,min(error_max2,3),'filled',DisplayName="Mixed")
scatter(Nt_list*dt,min(error_max3,3),'filled',DisplayName="Absorbing")
legend()
xlabel("$T$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Error analysis for various $T$ (cut-off at 3)","Interpreter","latex","FontSize",18)
hold off



%% Analysis of the time step dt

wh1 = discrete_laplace_transform(U1,t_list,10,"Trapezoidal");
[Ah1, submesh] = inverse_solver(wh1, mesh, 10);
a_exact = afun(submesh.vtx);

dt_list = dt:dt:100*dt;

error_quad1 = zeros(size(dt_list));
error_quad2 = zeros(size(dt_list));
error_quad3 = zeros(size(dt_list));

error_max1 = zeros(size(dt_list));
error_max2 = zeros(size(dt_list));
error_max3 = zeros(size(dt_list));


for k = 1:size(dt_list,2)

    wh1 = discrete_laplace_transform(U1(:,1:k:end),t_list(1:k:end),s,"Trapezoidal");
    [Ah1, ~] = inverse_solver(wh1, mesh, s);
    Ah1 = boundary_post_processing(Ah1,submesh);

    wh2 = discrete_laplace_transform(U2(:,1:k:end),t_list(1:k:end),s,"Trapezoidal");
    [Ah2, ~] = inverse_solver(wh2, mesh, s);
    Ah2 = boundary_post_processing(Ah2,submesh);

    wh3 = discrete_laplace_transform(U3(:,1:k:end),t_list(1:k:end),s,"Trapezoidal");
    [Ah3, ~] = inverse_solver(wh3, mesh, s);
    Ah3 = boundary_post_processing(Ah3,submesh);

    error_quad1(k) = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad2(k) = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_quad3(k) = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));

    error_max1(k) = max(abs(Ah1 - a_exact));
    error_max2(k) = max(abs(Ah2 - a_exact));
    error_max3(k) = max(abs(Ah3 - a_exact));
end


%% Plot the results

figure
scatter(dt_list,min(error_quad1,2),'filled',DisplayName="Neumann")
hold on
scatter(dt_list,min(error_quad2,2),'filled',DisplayName="Mixed")
scatter(dt_list,min(error_quad3,2),'filled',DisplayName="Absorbing")
legend()
xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18)
title("Error analysis for various $\mathrm{d}t$ (cut-off at 2)","Interpreter","latex","FontSize",18)
hold off

figure
scatter(dt_list,min(error_max1,3),'filled',DisplayName="Neumann")
hold on
scatter(dt_list,min(error_max2,3),'filled',DisplayName="Mixed")
scatter(dt_list,min(error_max3,3),'filled',DisplayName="Absorbing")
legend()
xlabel("$\mathrm{d}t$","Interpreter","latex","FontSize",18)
ylabel("$||a_h - a||_{\infty}$","Interpreter","latex","FontSize",18)
title("Error analysis for various $\mathrm{d}t$ (cut-off at 3)","Interpreter","latex","FontSize",18)
hold off
