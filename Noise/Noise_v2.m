%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the robustness to the noise using time-dependent data.


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

s = 5;                                                  % Pseudo-frequency

% Noise parameters

noise_lvl = 3/100;                                      % Noise level

alph = 2*rand-1;                                        % Random variable uniformly distributed


%% Generation of the unoised data

[U1, mesh, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);

[U2, ~, ~] = forward_solver_method2_2(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);

[U3, ~, ~] = forward_solver_method2_3(h,dt,Nt,afun,pfun,2*pi/om,delta,theta);


wh1 = discrete_laplace_transform(U1,t_list,s,"Trapezoidal");
wh2 = discrete_laplace_transform(U2,t_list,s,"Trapezoidal");
wh3 = discrete_laplace_transform(U3,t_list,s,"Trapezoidal");


[Ah_star1, submesh] = inverse_solver(wh1, mesh, s);
[Ah_star2, ~] = inverse_solver(wh2, mesh, s);
[Ah_star3, ~] = inverse_solver(wh3, mesh, s);


Ah_star1 = boundary_post_processing(Ah_star1,submesh);
Ah_star2 = boundary_post_processing(Ah_star2,submesh);
Ah_star3 = boundary_post_processing(Ah_star3,submesh);


Vh2 = fem(submesh, 'P1');

a_exact = afun(submesh.vtx);


%% Homogeneous noise

% Simulation of a reconstruction without smoothing

U_til1 = U1*(1 + noise_lvl*alph);
U_til2 = U2*(1 + noise_lvl*alph);
U_til3 = U3*(1 + noise_lvl*alph);


wh_til1 = discrete_laplace_transform(U_til1,t_list,s,"Trapezoidal");
wh_til2 = discrete_laplace_transform(U_til2,t_list,s,"Trapezoidal");
wh_til3 = discrete_laplace_transform(U_til3,t_list,s,"Trapezoidal");


[Ah1, ~] = inverse_solver(wh_til1, mesh, s);
[Ah2, ~] = inverse_solver(wh_til2, mesh, s);
[Ah3, ~] = inverse_solver(wh_til3, mesh, s);


Ah1 = boundary_post_processing(Ah1,submesh);
Ah2 = boundary_post_processing(Ah2,submesh);
Ah3 = boundary_post_processing(Ah3,submesh);


figure
graph(Vh2, Ah1)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ using noised data from Model 1",'Interpreter','latex','FontSize',18)

figure
graph(Vh2, Ah2)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ using noised data from Model 2",'Interpreter','latex','FontSize',18)


figure
graph(Vh2, Ah3)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ using noised data from Model 3",'Interpreter','latex','FontSize',18)



% Analysis of the effect of the homogeneous noise

alph_list = -1:0.01:1;

error_list1 = zeros(size(alph_list));
error_list2 = zeros(size(alph_list));
error_list3 = zeros(size(alph_list));

rel_error_list1 = zeros(size(alph_list));
rel_error_list2 = zeros(size(alph_list));
rel_error_list3 = zeros(size(alph_list));


for k = 1:size(alph_list,2)

    U_til1 = U1*(1+noise_lvl*alph_list(k));
    U_til2 = U2*(1+noise_lvl*alph_list(k));
    U_til3 = U3*(1+noise_lvl*alph_list(k));

    wh_til1 = discrete_laplace_transform(U_til1,t_list,s,"Trapezoidal");
    wh_til2 = discrete_laplace_transform(U_til2,t_list,s,"Trapezoidal");
    wh_til3 = discrete_laplace_transform(U_til3,t_list,s,"Trapezoidal");

    [Ah1, ~] = inverse_solver(wh_til1, mesh, s);
    [Ah2, ~] = inverse_solver(wh_til2, mesh, s);
    [Ah3, ~] = inverse_solver(wh_til3, mesh, s);

    Ah1 = boundary_post_processing(Ah1,submesh);
    Ah2 = boundary_post_processing(Ah2,submesh);
    Ah3 = boundary_post_processing(Ah3,submesh);

    error_list1(k) = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_list2(k) = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));
    error_list3(k) = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));

    rel_error_list1(k) = sqrt(sum((Ah1-Ah_star1).^2))/sqrt(sum(Ah_star1.^2));
    rel_error_list2(k) = sqrt(sum((Ah2-Ah_star2).^2))/sqrt(sum(Ah_star2.^2));
    rel_error_list3(k) = sqrt(sum((Ah3-Ah_star3).^2))/sqrt(sum(Ah_star3.^2));

end



figure
plot(alph_list, error_list1)
hold on
yline(error_list1(alph_list==0),'r--');
hold off
legend(["","Without noise"])
title("Analysis of the effect of the simple noise with $\delta =$"+sprintf('%1.0f', 100*noise_lvl)+" $\%$",'Interpreter','latex','FontSize',18);
xlabel("\alpha","FontSize",15);
ylabel("$\frac{||a_{h, noise} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);


figure
plot(alph_list, error_list2)
hold on
yline(error_list2(alph_list==0),'r--');
hold off
legend(["","Without noise"])
title("Analysis of the effect of the simple noise with $\delta =$"+sprintf('%1.0f', 100*noise_lvl)+" $\%$",'Interpreter','latex','FontSize',18);
xlabel("\alpha","FontSize",15);
ylabel("$\frac{||a_{h, noise} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);


figure
plot(alph_list, error_list3)
hold on
yline(error_list3(alph_list==0),'r--');
hold off
legend(["","Without noise"])
title("Analysis of the effect of the simple noise with $\delta =$"+sprintf('%1.0f', 100*noise_lvl)+" $\%$",'Interpreter','latex','FontSize',18);
xlabel("\alpha","FontSize",15);
ylabel("$\frac{||a_{h, noise} - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);


figure
plot(alph_list, rel_error_list1)
hold on
yline(rel_error_list1(alph_list==0),'r--');
hold off
title("Analysis of the effect of the simple noise with $\delta =$"+sprintf('%1.0f', 100*noise_lvl)+" $\%$",'Interpreter','latex','FontSize',18);
xlabel("\alpha","FontSize",15);
ylabel("$\frac{||a_{h, noise} - a_{h, unoised}||_2}{||a_{h, unoised}||_2}$","Interpreter","latex","FontSize",18);


figure
plot(alph_list, rel_error_list2)
hold on
yline(rel_error_list2(alph_list==0),'r--');
hold off
title("Analysis of the effect of the simple noise with $\delta =$"+sprintf('%1.0f', 100*noise_lvl)+" $\%$",'Interpreter','latex','FontSize',18);
xlabel("\alpha","FontSize",15);
ylabel("$\frac{||a_{h, noise} - a_{h, unoised}||_2}{||a_{h, unoised}||_2}$","Interpreter","latex","FontSize",18);


figure
plot(alph_list, rel_error_list3)
hold on
yline(rel_error_list3(alph_list==0),'r--');
hold off
title("Analysis of the effect of the simple noise with $\delta =$"+sprintf('%1.0f', 100*noise_lvl)+" $\%$",'Interpreter','latex','FontSize',18);
xlabel("\alpha","FontSize",15);
ylabel("$\frac{||a_{h, noise} - a_{h, unoised}||_2}{||a_{h, unoised}||_2}$","Interpreter","latex","FontSize",18);


%% Matrix noise

% Parameter for the method
U = U3;

Ah_star = Ah_star3;

error_0 = sqrt(sum((Ah_star-a_exact).^2))/sqrt(sum(a_exact.^2));

% Generation of noise
rng(1);
Alph = 2*rand(size(U))-1;

U_til = U .* (1+noise_lvl*Alph);

wh_til = discrete_laplace_transform(U_til,t_list,s,"Trapezoidal");


% Without smoothing
[Ah1, ~] = inverse_solver(wh_til, mesh, s);

error1 = sqrt(sum((Ah1-a_exact).^2))/sqrt(sum(a_exact.^2));

rel_error1 = sqrt(sum((Ah1-Ah_star).^2))/sqrt(sum(Ah_star.^2));


% 2D-pre-smoothing
[~,~,vh] = mesh2surface(wh_til,mesh);

% Smoothing method
vh_smooth = smoothdata2(vh,'gaussian');           % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'

v = reshape(vh_smooth,[],1);

[Ah2, ~] = inverse_solver(v, mesh, s);
Ah2 = boundary_post_processing(Ah2,submesh);

error2 = sqrt(sum((Ah2-a_exact).^2))/sqrt(sum(a_exact.^2));

rel_error2 = sqrt(sum((Ah2-Ah_star).^2))/sqrt(sum(Ah_star.^2));


% 2D-post-smoothing
[~,~,z] = mesh2surface(Ah2,submesh);

z = smoothdata2(z,'loess');                     % 'movmean' 'movmedian' 'gaussian' 'lowess' 'loess'

Ah3 = reshape(z, [], 1);

error3 = sqrt(sum((Ah3-a_exact).^2))/sqrt(sum(a_exact.^2));

rel_error3 = sqrt(sum((Ah3-Ah_star).^2))/sqrt(sum(Ah_star.^2));


% Smoothing at each time step
U_smooth = zeros(size(U_til));

for k = 1:size(U_til,2)
    [~,~,vh] = mesh2surface(U_til(:,k),mesh);

    vh_smooth = smoothdata2(vh,'gaussian');

    v = reshape(vh_smooth,[],1);

    U_smooth(:,k) = v;
end

wh_smooth = discrete_laplace_transform(U_smooth,t_list,s,"Trapezoidal");

[Ah4, ~] = inverse_solver(wh_smooth, mesh, s);
Ah4 = boundary_post_processing(Ah4,submesh);

[~,~,z] = mesh2surface(Ah4,submesh);
z = smoothdata2(z,'loess');
Ah5 = reshape(z, [], 1);

error4 = sqrt(sum((Ah4-a_exact).^2))/sqrt(sum(a_exact.^2));

rel_error4 = sqrt(sum((Ah4-Ah_star).^2))/sqrt(sum(Ah_star.^2));

error5 = sqrt(sum((Ah5-a_exact).^2))/sqrt(sum(a_exact.^2));

rel_error5 = sqrt(sum((Ah5-Ah_star).^2))/sqrt(sum(Ah_star.^2));


% Plot the results
figure
graph(Vh2, Ah_star)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ without noise",'Interpreter','latex','FontSize',18)

figure
graph(Vh2, Ah1)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ without smoothing",'Interpreter','latex','FontSize',18)

figure
graph(Vh2, Ah2)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ with 2D-pre-smoothing",'Interpreter','latex','FontSize',18)

figure
graph(Vh2, Ah3)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ with 2D-pre/post-smoothing",'Interpreter','latex','FontSize',18)

figure
graph(Vh2, Ah4)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ with 2D/time-pre-smoothing",'Interpreter','latex','FontSize',18)

figure
graph(Vh2, Ah5)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ with 2D/time-pre/post-smoothing",'Interpreter','latex','FontSize',18)
