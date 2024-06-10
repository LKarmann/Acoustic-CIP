%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Solves the inverse problem without the knowledge of s but assuming the
% boundary condition is known.


%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Fonctions\;


%% Definition of parameters

h = 2^(-5);                                             % Mesh size

s_star = 10.2;                                          % True pseudo-frequency

s_inf = 0.5;

s_sup = 30;

s_step = 0.5;

tol = 0.001;

om = 80;                                                % Frequency of the boudary source

c = (1-exp(-2*pi*s_star/om))/((1+(s_star/om)^2)*om);    % Boundary source

afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);



%% Generation of the measurement

[wh, mesh] = forward_solver_method1(h, s_star, afun, pfun);



%% Solving the blind inverse problem

s_list = s_inf:s_step:s_sup;

error_list = zeros(size(s_list));

true_error_list = zeros(size(s_list));



function Y = Anfun(X, Anh, submesh)

Y = zeros([size(X,1), 1]);

for j = 1:size(X,1)
    if min(X(j,1), X(j,2)) > 0 && max(X(j,1), X(j,2)) < 1

        [~,I] = min((submesh.vtx(:,1)-X(j,1)).^2 + (submesh.vtx(:,2)-X(j,2)).^2);

        Y(j) = Anh(I);

    else

        Y(j) = 1;

    end
end

end

k = 1;


while k <= size(s_list,2)

    s = s_list(k);

    [Anh, submesh] = inverse_solver(wh, mesh, s);


    % Calculation of the true error
    A = afun(submesh.vtx);

    true_error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));


    % Calculation of the relative error
    wnh = forward_solver_method1(h,s,@(X) Anfun(X,Anh,submesh),pfun);

    error_list(k) = sqrt(sum((wnh-wh).^2))/sqrt(sum(wh.^2));

    k = k+1;

end




% The first proposal (Inefficient at all)
% % Initialisation k=1 and k=2
% 
% s = s_list(1);
% 
% [Anh, submesh] = inverse_solver(wh, mesh, s);
% 
% A = afun(submesh.vtx);
% 
% true_error_list(1) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
% 
% Anm1h = Anh;
% 
% 
% 
% s = s_list(2);
% 
% [Anh, submesh] = inverse_solver(wh, mesh, s);
% 
% A = afun(submesh.vtx);
% 
% true_error_list(2) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
% 
% error_list(1) = 2*sqrt(sum((Anh-Anm1h).^2))/sqrt(sum(Anh.^2));
% 
% error_list(2) = sqrt(sum((Anh-Anm1h).^2))/sqrt(sum(Anh.^2));
% 
% Anm1h = Anh;
% 
% 
% k = 2;
% 
% 
% 
% while k <= size(s_list,2)-1 %&& error_list(k) <= error_list(k-1) && error_list(k) > tol
% 
%     k = k+1;
% 
%     s = s_list(k);
% 
%     [Anh, submesh] = inverse_solver(wh, mesh, s);
% 
%     A = afun(submesh.vtx);
% 
%     true_error_list(k) = sqrt(sum((Anh-A).^2))/sqrt(sum(A.^2));
% 
%     error_list(k) = sqrt(sum((Anh-Anm1h).^2))/sqrt(sum(Anh.^2));
% 
%     Anm1h = Anh;
% 
% end



%% Plot the results

figure
plot(s_list, error_list)
title("Relative quadratic error on $w_h$","Interpreter","latex");
xlabel("s","Interpreter","latex");
ylabel("$\frac{|| w_{n,h} - w_h ||_2}{||w_h||_2}$","Interpreter","latex");
axis([0,30,0,1]);

figure
plot(s_list, true_error_list)
title("Relative quadratic error on $a$","Interpreter","latex");
xlabel("s","Interpreter","latex");
ylabel("$\frac{|| a_n - a ||_2}{||a||_2}$","Interpreter","latex");
axis([0,30,0,1]);




% %% Solving using blind_inverse_solver
% 
% [Anh, submesh, s, error_k] = blind_inverse_solver(wh, mesh, h, s_inf, s_sup, s_step, tol, pfun, "increase");
% 
% Vh2 = fem(submesh, 'P1');
% 
% figure
% graph(Vh2, Anh)
% title('Reconstruction of $a_{approx}$','interpreter','latex');
% axis([-0.1, 1.1, -0.1, 1.1, 0, 4]);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% zlabel('$z$','interpreter','latex');