%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Solves the Acoustic Coefficient Inverse Problem with Algorithm 2 and
% Method 1.


%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions\;


%% Definition of parameters

h = 2^(-5);                                             % Mesh size

s = 10;                                                 % Pseudo-frequency

om = 80;                                                % Frequency of the boudary source

c = (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om);              % Boundary source

Index = 1;                                              % Index for II

Extraction = 1;                                         % Extraction for inverse_solver_2

Proposal = 3;                                           % Proposal value

% afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);



%% Solve the inverse problem through other functions

[wh, mesh] = forward_solver_method1(h, s, afun, pfun);


if Index == 1

    [Anh, submesh] = inverse_solver_1(wh, mesh, h, s, Proposal);

else

    [Anh, submesh] = inverse_solver_2(wh, mesh, h, s, pfun, Extraction);

end


% Direct resolution
% Wnh = wh(mesh.vtx(:,1) >=0 & mesh.vtx(:,2) >=0 & mesh.vtx(:,1) <=1 & mesh.vtx(:,2) <=1);
% 
% [x,y,z]=submeshsurface(Wnh,submesh);
% 
% delz = 4*del2(z,x,y);
% 
% 
% delw = reshape(delz,[],1);
% 
% Anh = delw ./ Wnh /s^2;


%% Calculation of the exact solution on the same mesh and the error

Vh2 = fem(submesh, 'P1');

A = afun(submesh.vtx);

error_quad = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));

error_max = max(abs(Anh - A));

B = unique(abs(Anh - A));
thres = 0.01;

error_max_thres = B(end-round(thres*size(B,1)));



% Plot the results

figure
graph(Vh2, A)
title('Input value of $a$','interpreter','latex');
axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');

figure
graph(Vh2, Anh)
title('Reconstruction of $a_{approx}$','interpreter','latex');
axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');

figure
graph(Vh2, Anh - A)
title('Difference between $a$ and $a_{approx}$','interpreter','latex');
axis([-0.1, 1.1, -0.1, 1.1, -1, 2]);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$','interpreter','latex');