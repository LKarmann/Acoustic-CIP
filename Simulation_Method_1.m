%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Simulations and figures of the results from Forward solver and
% reconstruction of the acoustic coefficient using Method 1 and Inverse
% solvers version 1 and 2 (see Appendix).


%% Preamble

clear;
close all;
format long;

addpath \gypsilab-master\gypsilab-master\openMsh;
addpath \gypsilab-master\gypsilab-master\openDom;
addpath \gypsilab-master\gypsilab-master\openFem;
addpath \Functions;


%% Definition of the parameters

h = 2^(-5);                                             % Mesh size

s = 5;                                                  % Pseudo-frequency

om = 80;                                                % Frequency of the boudary source

c = (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om);              % Laplace transform of p(t) at s

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/3);              % Boundary source

% Test 1: 1 inclusion
afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% Test 2: 2 inclusions
% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

Index = 1;                                              % Selection of the inverse solver between 1 and 2

Extraction = 1;                                         % Extraction parameter for inverse_solver_2

Proposal = 3;                                           % Proposal for the calculation of the normal derivative


%% Solve the forward problem using Method 1

[wh, mesh] = forward_solver_method1(h,s,afun,pfun);


%% Solve the inverse problem

if Index == 1
    [Ah, submesh] = inverse_solver_1(wh,mesh,h,s,Proposal);
else
    [Ah, submesh] = inverse_solver_2(wh,mesh,h,s,pfun,Extraction);
end

%% Post-processing the result

Ah_2 = boundary_post_processing(Ah,submesh);


%% Plot the result

Vh = fem(submesh,'P1');

a_exact = afun(submesh.vtx);

figure
graph(Vh, a_exact)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a(x,y)$ for Test 1",'Interpreter','latex','FontSize',18)

figure
graph(Vh, Ah)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ for Test 1 at $s=5$",'Interpreter','latex','FontSize',18)

figure
graph(Vh, Ah_2)
colorbar()
xlabel("$x$",'Interpreter','latex','FontSize',18)
ylabel("$y$",'Interpreter','latex','FontSize',18)
title("$a_h(x,y)$ for Test 1 at $s=5$",'Interpreter','latex','FontSize',18)


%% Calculation of the relative quadratic error and of the maximal error

error_quad = sqrt(sum((Ah-a_exact).^2))/sqrt(sum((a_exact).^2));

error_max = max(abs(Ah-a_exact));

error_quad_2 = sqrt(sum((Ah_2-a_exact).^2))/sqrt(sum((a_exact).^2));

error_max_2 = max(abs(Ah_2-a_exact));
