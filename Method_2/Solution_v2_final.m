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

dt = 0.002;                                             % Time step (in seconds)

Nt = 1000;                                              % Number of time steps

t_list = 0:dt:Nt*dt;                                    % Time discretization

om = 80;                                                % Frequency of the boudary source

pfun = @(t) sin(om*t) * (om * t <= 2*pi);               % Boundary source

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

delta = 1;                                              % Parameter of the numerical integration

theta = 1/2;                                            % Parameter of the numerical integration



% Inverse problem
s = 10;                                                 % Pseudo-frequency

pfun1 = @(X) (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om) * (abs(X(:,2) - 1.5) < h/3);

%% Solving the forward problem

[wh1, mesh1] = forward_solver_method1(h,s,afun,pfun1);

[U_list2, mesh2, t_list2] = forward_solver_method2_1(h,dt,Nt,afun,pfun,delta,theta);

[U_list3, mesh3, t_list3] = forward_solver_method2_2(h,dt,Nt,afun,pfun,om,delta,theta,1);

[U_list4, mesh4, t_list4] = forward_solver_method2_3(h,dt,Nt,afun,pfun,om,delta,theta,1);



% %% Movie of the solution
% 
% Vh = fem(mesh3,'P1');
% 
% figure
% hold on
% pause(3)
% for k = 1:5:size(U_list3,2)
%     clf
%     graph(Vh,U_list3(:,k))
%     drawnow
% end
% hold off




%% Calculation of the Laplace transform

wh2 = discrete_laplace_transform(U_list2,t_list2,s,"LeftRectangle");

wh3 = discrete_laplace_transform(U_list3,t_list3,s,"LeftRectangle");

wh4 = discrete_laplace_transform(U_list4,t_list4,s,"LeftRectangle");




%% Solving the inverse problem

[Anh1, submesh] = inverse_solver(wh1,mesh1,s);
[Anh2, ~] = inverse_solver(wh2,mesh2,s);
[Anh3, ~] = inverse_solver(wh3,mesh3,s);
[Anh4, ~] = inverse_solver(wh4,mesh4,s);

Anh1_2 = Anh1.*(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2))...
    +1-(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2));
Anh2_2 = Anh2.*(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2))...
    +1-(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2));
Anh3_2 = Anh3.*(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2))...
    +1-(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2));
Anh4_2 = Anh4.*(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2))...
    +1-(submesh.vtx(:,1)>=h/sqrt(2)&(1-submesh.vtx(:,1))>=h/sqrt(2)&submesh.vtx(:,2)>=h/sqrt(2)&(1-submesh.vtx(:,2))>=h/sqrt(2));

Vh2 = fem(submesh,'P1');

A = afun(submesh.vtx);



% %% Plot the Results
% 
% figure
% graph(Vh2, A)
% title('Input value of $a$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2, Anh5)
% title('Reconstruction of $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');
% 
% figure
% graph(Vh2, Anh5 - A)
% title('Difference between $a$ and $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, -1, 2]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');



%% Calculation of the relative error

error_quad1 = sqrt(sum((Anh1-A).^2))/sqrt(sum(A.^2));

B = unique(abs(Anh1 - A));
thres = 0.01;

error_max1 = B(end);

error_max_thres1 = B(end-round(thres*size(B,1)));



error_quad2 = sqrt(sum((Anh2-A).^2))/sqrt(sum(A.^2));

B = unique(abs(Anh2 - A));
thres = 0.01;

error_max2 = B(end);

error_max_thres2 = B(end-round(thres*size(B,1)));



error_quad3 = sqrt(sum((Anh3-A).^2))/sqrt(sum(A.^2));

B = unique(abs(Anh3 - A));
thres = 0.01;

error_max3 = B(end);

error_max_thres3 = B(end-round(thres*size(B,1)));



error_quad4 = sqrt(sum((Anh4-A).^2))/sqrt(sum(A.^2));

B = unique(abs(Anh4 - A));
thres = 0.01;

error_max4 = B(end);

error_max_thres4 = B(end-round(thres*size(B,1)));



%% Post-processing

error_quad1_2 = sqrt(sum((Anh1_2-A).^2))/sqrt(sum(A.^2));

error_max1_2 = max(abs(Anh1_2-A));



error_quad2_2 = sqrt(sum((Anh2_2-A).^2))/sqrt(sum(A.^2));

error_max2_2 = max(abs(Anh2_2-A));



error_quad3_2 = sqrt(sum((Anh3_2-A).^2))/sqrt(sum(A.^2));

error_max3_2 = max(abs(Anh3_2-A));



error_quad4_2 = sqrt(sum((Anh4_2-A).^2))/sqrt(sum(A.^2));

error_max4_2 = max(abs(Anh4_2-A));
