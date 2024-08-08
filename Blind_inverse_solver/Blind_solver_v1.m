%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the principle of the blind inverse solver and analysis of
% the error on the boundary condition (see Appendix).


%% Preamble

clear;
close all;
format long;

addpath \gypsilab-master\gypsilab-master\openMsh;
addpath \gypsilab-master\gypsilab-master\openDom;
addpath \gypsilab-master\gypsilab-master\openFem;
addpath Functions;


%% Definition of the parameters

h = 2^(-5);                                                 % Mesh size

s_star = 10.2;                                              % Pseudo-frequency of reference

om = 80;                                                    % Frequency of the boudary source

c_star = (1-exp(-2*pi*s_star/om))/((1+(s_star/om)^2)*om);   % Laplace transform of p(t) at s_star

pfun = @(X) c_star * (abs(X(:,2) - 1.5) < h/2);             % Boundary source

s_inf = 5;                                                  % Inf value of s

s_sup = 20;                                                 % Sup value of s

s_step = 0.1;                                               % Step size for s

tol = 0.001;                                                % Tolerance for the blind inverse solver

noise_lvl = 10/100;                                         % Noise level

% Test 1: 1 inclusion
afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% Test 2: 2 inclusions
% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);


%% Generation of the measurement and of the reference reconstruction

[wh_star, mesh] = forward_solver_method1(h, s_star, afun, pfun);

[Ah_star, submesh] = inverse_solver(wh_star,mesh,s_star);
Ah_star = boundary_post_processing(Ah_star,submesh);

a_exact = afun(submesh.vtx);

error_ref = sqrt(sum((Ah_star-a_exact).^2))/sqrt(sum(a_exact.^2));


%% Analysis of the error for various s

s_list = s_inf:s_step:s_sup;

error_list = zeros(size(s_list));

error_a_list = zeros(size(s_list));


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

    [Ah, submesh] = inverse_solver(wh_star, mesh, s);
    Ah = boundary_post_processing(Ah,submesh);

    % Calculation of the error on a
    error_a_list(k) = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));

    % Calculation of the error on w
    wh = forward_solver_method1(h,s,@(X) Afun(X,Ah,submesh),pfun);

    error_list(k) = sqrt(sum((wh-wh_star).^2))/sqrt(sum(wh_star.^2));
end


%% Plot the results

figure
plot(s_list, error_list)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| w_k - w(s_*) ||_2}{||w(s_*)||_2}$","Interpreter","latex","FontSize",18);
axis([0,30,0,1]);

figure
plot(s_list, error_a_list)
hold on
xline(s_star,'-r',"$s^*$","Interpreter","latex",'FontSize',14)
hold off
title("Reconstruction without knowledge of $s$","Interpreter","latex","FontSize",18);
xlabel("$s$","Interpreter","latex","FontSize",18);
ylabel("$\frac{|| a_k - a ||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
axis([0,30,0,1]);


%% Solve using blind_inverse_solver

[Ah, submesh, s, error_k] = blind_inverse_solver(wh_star,mesh,h,s_inf,s_sup,s_step,tol,pfun,"increase");
Ah = boundary_post_processing(Ah,submesh);

Vh2 = fem(submesh,'P1');

figure
graph(Vh2, Ah)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ using the blind inverse solver",'Interpreter','latex','FontSize',18)


%% Blind inverse solver with an error on c_star

c_list = linspace(c_star*(1-noise_lvl), c_star*(1+noise_lvl), 15);

s_list = zeros(size(c_list));
error_list = zeros(size(c_list));


for k = 1:size(c_list,2)
    pfun = @(X) c_list(k) * (abs(X(:,2) - 1.5) < h/2);

    [Ah, submesh, s, error] = blind_inverse_solver(wh_star,mesh,h,s_inf,s_sup,s_step,tol,pfun,"increase");
    Ah = boundary_post_processing(Ah,submesh);

    s_list(k) = s;

    error_list(k) = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));

end

% Simple linear regression
alpha_0 = (s_list-s_star)*(c_list-c_star)'/sum((c_list-c_star).^2)*c_star/s_star;

figure
scatter(c_list,s_list,'blue','filled')
hold on
plot(c_list,s_star*c_list/c_star)
plot(c_list,alpha_0*s_star*(c_list-c_star)/c_star + s_star)
hold off
legend(["$s_{reconstructed}$","$s = \frac{s_*}{\mathcal{L}(p)_*} \mathcal{L}(p) $",...
    "$s = \alpha \frac{s_*}{\mathcal{L}(p)_*} (\mathcal{L}(p) - \mathcal{L}(p)_*) + s_* $"],'Interpreter','latex','Location','northwest')
title("Reconstruction of $s$ with noise on $\mathcal{L}(p)_*$","Interpreter","latex","FontSize",18);
text(0.0072,9.2,"$\alpha = $"+sprintf('%1.6f',alpha_0),'Interpreter','latex')
xlabel("$\mathcal{L}(p)$",'Interpreter','latex','FontSize',18);
ylabel("$s$",'Interpreter','latex','FontSize',18);


figure
plot(c_list, error_list)
title("Blind inverse solver with noise on $\mathcal{L}(p)_*$","Interpreter","latex","FontSize",18);
xline(c_star,'r',"$\mathcal{L}(p)_*$",'Interpreter','latex');
xlabel("$\mathcal{L}(p)$",'Interpreter','latex','FontSize',18);
ylabel("$\frac{||a_h - a||_2}{||a||_2}$",'Interpreter','latex','FontSize',18);


%% Blind inverse solver using a noised value of s_star to calculate c_star

s_noise_list = linspace(s_star*(1-noise_lvl), s_star*(1+noise_lvl), 15);

s_list = zeros(size(s_noise_list));
error_list = zeros(size(s_noise_list));


for k = 1:size(s_noise_list,2)
    c = (1-exp(-2*pi*s_noise_list(k)/om))/((1+(s_noise_list(k)/om)^2)*om);

    pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);

    [Ah, submesh, s, error] = blind_inverse_solver(wh_star,mesh,h,s_inf,s_sup,s_step,tol,pfun,"increase");
    Ah = boundary_post_processing(Ah,submesh);

    s_list(k) = s;

    error_list(k) = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));

end

% Simple linear regression
alpha_0 = (s_list-s_star)*(s_noise_list-s_star)'/sum((s_noise_list-s_star).^2);

figure
scatter(s_noise_list,s_list,'blue','filled')
hold on
plot(s_noise_list,s_noise_list)
plot(s_noise_list,alpha_0*(s_noise_list-s_star) + s_star)
hold off
legend(["","Identity",...
    "$s = \alpha (s - s_*) + s_* $"],'Interpreter','latex','Location','northwest')
title("Reconstruction of $s$ for various $\mathcal{L}(p)(s)$","Interpreter","latex","FontSize",18);
text(11,9.2,"$\alpha = $"+sprintf('%1.6f',alpha_0),'Interpreter','latex')
xlabel("$s$",'Interpreter','latex','FontSize',18);
ylabel("$s_{reconstructed}$",'Interpreter','latex','FontSize',18);


figure
plot(s_noise_list, error_list)
title("Blind inverse solver for various $\mathcal{L}(p)(s)$","Interpreter","latex","FontSize",18);
xline(s_star,'r',"$s_*$",'Interpreter','latex');
xlabel("$s$",'Interpreter','latex','FontSize',18);
ylabel("$\frac{||a_h - a||_2}{||a||_2}$",'Interpreter','latex','FontSize',18);
