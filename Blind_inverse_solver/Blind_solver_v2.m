%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the principle of the blind inverse solver using time-
% dependent solvers.


%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions;


%% Definition of the parameters

h = 2^(-5);                                                 % Mesh size

dt = 0.001;                                                 % Time step (in seconds)

Nt = 10000;                                                 % Number of time steps

t_list = 0:dt:Nt*dt;                                        % Time discretization

om = 80;                                                    % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);                   % Boundary source

t_1 = 2*pi/om;                                              % Duration of the external source

s_star = 10.2;                                              % Pseudo-frequency of reference

c_star = (1-exp(-2*pi*s_star/om))/((1+(s_star/om)^2)*om);   % Laplace transform of p(t) at s_star

pfunbis = @(X) c_star * (abs(X(:,2) - 1.5) < h/2);          % Boundary source

s_inf = 5;                                                  % Inf value of s

s_sup = 20;                                                 % Sup value of s

s_step = 0.5;                                               % Step size for s

tol = 0.001;                                                % Tolerance for the blind inverse solver

noise_lvl = 10/100;                                         % Noise level

% Test 1: 1 inclusion
% afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% Test 2: 2 inclusions
afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);


%% Generation of the measurement and of the reference reconstruction

[Uh_star1, mesh, t_list] = forward_solver_method2_1(h, dt, Nt, afun, pfun, 1/2, 1/2);
[Uh_star2, ~, ~] = forward_solver_method2_2(h, dt, Nt, afun, pfun, t_1, 1/2, 1/2);
[Uh_star3, ~, ~] = forward_solver_method2_3(h, dt, Nt, afun, pfun, t_1, 1/2, 1/2);

wh_star1 = discrete_laplace_transform(Uh_star1,t_list,s_star,"Trapezoidal");
wh_star2 = discrete_laplace_transform(Uh_star2,t_list,s_star,"Trapezoidal");
wh_star3 = discrete_laplace_transform(Uh_star3,t_list,s_star,"Trapezoidal");

[Ah_star1, submesh] = inverse_solver(wh_star1,mesh,s_star);
[Ah_star2, ~] = inverse_solver(wh_star2,mesh,s_star);
[Ah_star3, ~] = inverse_solver(wh_star3,mesh,s_star);

Ah_star1 = boundary_post_processing(Ah_star1,submesh);
Ah_star2 = boundary_post_processing(Ah_star2,submesh);
Ah_star3 = boundary_post_processing(Ah_star3,submesh);

a_exact = afun(submesh.vtx);

error_ref1 = sqrt(sum((Ah_star1-a_exact).^2))/sqrt(sum(a_exact.^2));
error_ref2 = sqrt(sum((Ah_star2-a_exact).^2))/sqrt(sum(a_exact.^2));
error_ref3 = sqrt(sum((Ah_star3-a_exact).^2))/sqrt(sum(a_exact.^2));


%% Analysis of the error for various s

s_list = s_inf:s_step:s_sup;

error_list1 = zeros(size(s_list));
error_list2 = zeros(size(s_list));
error_list3 = zeros(size(s_list));

error_a_list1 = zeros(size(s_list));
error_a_list2 = zeros(size(s_list));
error_a_list3 = zeros(size(s_list));


% Auxiliar function
function Y = Afun(X, Ah, submesh)
Y = zeros([size(X,1), 1]);
for j = 1:size(X,1)
    if min(X(j,1), X(j,2)) > 0 && max(X(j,1), X(j,2)) < 1
        [~,I] = min((submesh.vtx(:,1)-X(j,1)).^2 + (submesh.vtx(:,2)-X(j,2)).^2);
        Y(j) = Ah(I);
    else
        Y(j) = 1;
    end
end
end


for k = 1:size(s_list,2)

    s = s_list(k);

    [Ah1, submesh] = inverse_solver(wh_star1, mesh, s);
    [Ah2, ~] = inverse_solver(wh_star2, mesh, s);
    [Ah3, ~] = inverse_solver(wh_star3, mesh, s);

    Ah1 = boundary_post_processing(Ah1,submesh);
    Ah2 = boundary_post_processing(Ah2,submesh);
    Ah3 = boundary_post_processing(Ah3,submesh);

    % Calculation of the error on a
    error_a_list1(k) = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_a_list2(k) = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_a_list3(k) = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));

    % Calculation of the error on w
    [Uh1, mesh, t_list] = forward_solver_method2_1(h,dt,Nt,@(X) Afun(X,Ah1,submesh),pfun,1/2,1/2);
    [Uh2, ~, ~] = forward_solver_method2_2(h,dt,Nt,@(X) Afun(X,Ah2,submesh),pfun,t_1,1/2,1/2);
    [Uh3, ~, ~] = forward_solver_method2_3(h,dt,Nt,@(X) Afun(X,Ah3,submesh),pfun,t_1,1/2,1/2);
    
    wh1 = discrete_laplace_transform(Uh1,t_list,s,"Trapezoidal");
    wh2 = discrete_laplace_transform(Uh2,t_list,s,"Trapezoidal");
    wh3 = discrete_laplace_transform(Uh3,t_list,s,"Trapezoidal");

    error_list1(k) = sqrt(sum((wh1-wh_star1).^2))/sqrt(sum(wh_star1.^2));
    error_list2(k) = sqrt(sum((wh2-wh_star2).^2))/sqrt(sum(wh_star2.^2));
    error_list3(k) = sqrt(sum((wh3-wh_star3).^2))/sqrt(sum(wh_star3.^2));
end


%% Plot the results

figure
plot(s_list, error_list1)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| w_k - w(s_*) ||_2}{||w(s_*)||_2}$","Interpreter","latex","FontSize",18);
axis([5,20,0,1]);

figure
plot(s_list, error_a_list1)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| a_k - a ||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
axis([0,30,0,1]);

figure
plot(s_list, error_list2)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| w_k - w(s_*) ||_2}{||w(s_*)||_2}$","Interpreter","latex","FontSize",18);
axis([0,30,0,1]);

figure
plot(s_list, error_a_list2)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| a_k - a ||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
axis([0,30,0,1]);

figure
plot(s_list, error_list3)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| w_k - w(s_*) ||_2}{||w(s_*)||_2}$","Interpreter","latex","FontSize",18);
axis([0,30,0,1]);

figure
plot(s_list, error_a_list3)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| a_k - a ||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
axis([0,30,0,1]);


figure
plot(s_list, error_a_list1,DisplayName="Neumann");
hold on
plot(s_list, error_a_list2,DisplayName="Mixed");
plot(s_list, error_a_list3,DisplayName="Absorbing");
scatter(s_star,error_ref1,"filled",DisplayName="s_*");
title("Reconstruction with an unknown $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{||a_k - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
hold off

figure
plot(s_list, error_a_list1,DisplayName="Neumann");
hold on
plot(s_list, error_a_list2,DisplayName="Mixed");
plot(s_list, error_a_list3,DisplayName="Absorbing");
scatter(s_star,error_ref1,"filled",DisplayName="s_*");
axis([s_star*(1-noise_lvl),s_star*(1+noise_lvl),0,1]);
xticks([s_star*(1-noise_lvl), 9.5,10,10.5,11,s_star*(1+noise_lvl)])
xticklabels(["s_* - 10%", "9.5","10","10.5","11","s_*+10%"])
legend("Neumann","Mixed","Absorbing","s_*")
title("Reconstruction with an unknown $s$ (zoom-in)","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{||a_k - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
hold off


%% Solve using blind_inverse_solver_time

[Ah, submesh, s, error_k] = blind_inverse_solver_time(wh_star3,mesh,h,dt,Nt,s_inf,s_sup,s_step,tol,pfun,"Absorbing",t_1,"increase");
Ah = boundary_post_processing(Ah,submesh);

Vh2 = fem(submesh,'P1');

figure
graph(Vh2, Ah)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ using the blind inverse solver",'Interpreter','latex','FontSize',18)

error_quad = sqrt(sum((Ah-a_exact).^2))/sqrt(sum((a_exact).^2));

error_max = max(abs(Ah-a_exact));


%% Solve using blind_inverse_solver

[Ah, submesh, s, error_k] = blind_inverse_solver(wh_star3,mesh,h,s_inf,s_sup,s_step,tol,pfunbis,"increase");
Ah = boundary_post_processing(Ah,submesh);

Vh2 = fem(submesh,'P1');

figure
graph(Vh2, Ah)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ using the blind inverse solver",'Interpreter','latex','FontSize',18)