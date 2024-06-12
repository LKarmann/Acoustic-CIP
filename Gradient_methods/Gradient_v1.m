%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Simulations of the gradient method

%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions;


%% Definition of a and parameters

% Forward problem

h = 2^(-5);                                             % Mesh size

dt = 0.001;                                             % Time step (in seconds)

Nt = 2000;                                              % Number of time steps

om = 80;                                                % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

N = 7;

gamma_0 = 0.01;

q = 1;

tol = 0.95;

method = "linear";



%% Generation of the measurement

[u_list1, mesh_0, t_list] = forward_solver_method2_1(h,dt,Nt,afun,pfun,1/2,1/2);

step = mesh_0.stp;
dx = step(1);

bound_mask = abs(mesh_0.vtx(:,2)-1.5) <= dx/3;

% bound_mask = abs(mesh.vtx(:,1)+0.5) <= dx/3 | abs(mesh.vtx(:,2)+0.5) <= dx/3 | abs(mesh.vtx(:,1)-1.5) <= dx/3 | abs(mesh.vtx(:,2)-1.5) <= dx/3;

intern_mask = mesh_0.vtx(:,1) >= 0 & mesh_0.vtx(:,1) <= 1 & mesh_0.vtx(:,2) >= 0 & mesh_0.vtx(:,2) <= 1;

u_til = u_list1(bound_mask,:);

a_exact = afun(mesh_0.vtx);

a_0 = ones(size(a_exact));





%% Gradient Method

a_list_1 = gradient_1_1(u_til,a_0,1,0.1,N,mesh_0,dt,Nt,pfun,0,"None",1);

a_list_2 = gradient_1_1(u_til,a_0,0.01,1,N,mesh_0,dt,Nt,pfun,0,"Cutoff",1);

a_list_3 = gradient_1_1(u_til,a_0,0.01,1,N,mesh_0,dt,Nt,pfun,0,"Putup",1);



%% Adaptive Gradient Method


[u_list1, mesh, t_list] = forward_solver_method2_1(2^(-2),dt,Nt,afun,pfun,1/2,1/2);

step = mesh.stp;
dx = step(1);

bound_mask = abs(mesh.vtx(:,2)-1.5) <= dx/3;

% bound_mask = abs(mesh.vtx(:,1)+0.5) <= dx/3 | abs(mesh.vtx(:,2)+0.5) <= dx/3 | abs(mesh.vtx(:,1)-1.5) <= dx/3 | abs(mesh.vtx(:,2)-1.5) <= dx/3;

u_til = u_list1(bound_mask,:);

a_exact = afun(mesh.vtx);

a_0 = ones(size(a_exact));



a_list_4 = adaptive_gradient_1_1(u_til,a_0,1,0.1,0.95,N,mesh,dt,Nt,pfun,"linear",0,"Putup",1);

a_list_5 = adaptive_gradient_1_1(u_til,a_0,0.01,1,0.95,N,mesh,dt,Nt,pfun,"linear",0,"Cutoff",1);

a_list_6 = adaptive_gradient_1_1(u_til,a_0,0.01,1,0.95,N,mesh,dt,Nt,pfun,"linear",0,"Putup",1);



%% Ploting the results


Vh = fem(mesh, 'P1');


for k = 1:N+1

    figure
    graph_mesh(a_list_1(:,k),mesh_0,0)

end



for k = 1:N+1

    figure
    graph_mesh(a_list_5{k,1},a_list_5{k,2},1)
    axis([-0.5,1.5,-0.5,1.5,min(a_list_5{k,1})-0.1,max(a_list_5{k,1})+0.1])

end