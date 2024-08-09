%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Generation of the data used in Simulation 1 (see Appendix).


%% Preamble

clear;
close all;
format long;

addpath \gypsilab-master\gypsilab-master\openMsh;
addpath \gypsilab-master\gypsilab-master\openDom;
addpath \gypsilab-master\gypsilab-master\openFem;
addpath \Functions;


%% Definition of the general parameters

om = 80;                                                    % Frequency of the boudary source

% Test 1: 1 inclusion
afun = @(X) 1 + 2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

% Test 2: 2 inclusions
% afun = @(X) 1 + (2*exp(-((X(:,1)-0.5).^2 + (X(:,2)-0.7).^2)/0.001) + 3*exp(-((X(:,1)-0.2).^2 + (X(:,2)-0.6).^2)/0.001)) .* (min(X(:,1), X(:,2)) > 0 & max(X(:,1), X(:,2)) < 1);

thres = 0.01;                                               % Threshold when boundary_post_processing is not applied


%% Calculation for various parameters

C = {'Index', 'Proposal', 'h', 's', 'Ah', 'submesh', 'error_quad', 'error_max', 'error_max_thres'};

s_list = 0.5:0.5:30;

ind = 1;


% Index 1 = Resolution in Omega_1

Index = 1;                                                  % Selection of the inverse solver between 1 and 2

for k_h = 2:6
for k_s = 1:size(s_list,2)
for Proposal = 1:3

% Definition of parameters

ind = ind + 1;                                              % Index of the simulation

h = 2^(-k_h);                                               % Mesh size

s = s_list(k_s);                                            % Pseudo-frequency

c = (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om);                  % Laplace transform of p(t) at s

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);                  % Boundary source


% Reconstruction of the test function

[wh, mesh] = forward_solver_method1(h, s, afun, pfun);

[Ah, submesh] = inverse_solver_1(wh, mesh, h, s, Proposal);


% Computation of the exact solution on the same mesh and the error

a_exact = afun(submesh.vtx);

error_quad = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));

B = unique(abs(Ah - a_exact));

error_max = B(end);

error_max_thres = B(end-round(thres*size(B,1)));


C(ind,:) = {Index, Proposal, h, s, Ah, submesh, error_quad, error_max, error_max_thres};

end
end
end


% Index 2 = Resolution in Omega

Index = 2;

for k_h = 2:6
for k_s = 1:size(s_list,2)
for Extraction = 0:1

% Definition of parameters

ind = ind + 1;                                              % Index of the simulation

h = 2^(-k_h);                                               % Mesh size

s = s_list(k_s);                                            % Pseudo-frequency

c = (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om);                  % Laplace transform of p(t) at s

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);                  % Boundary source


% Solve the inverse problem through other functions

[wh, mesh] = forward_solver_method1(h, s, afun, pfun);

[Ah, submesh] = inverse_solver_2(wh, mesh, h, s, pfun, Extraction);


% Calculation of the exact solution on the same mesh and the error

a_exact = afun(submesh.vtx);

error_quad = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));

B = unique(abs(Ah - a_exact));

error_max = B(end);

error_max_thres = B(end-round(thres*size(B,1)));


C(ind,:) = {Index, Extraction, h, s, Ah, submesh, error_quad, error_max, error_max_thres};

end
end
end

% Save the data in a file

save Data_simulation_1.mat 'C'




%% Calculation for various parameters with boundary post-processing

C = {'Index', 'Proposal', 'h', 's', 'Ah', 'submesh', 'error_quad', 'error_max'};

s_list = 0.5:0.5:30;

ind = 1;


% Index 1 = Resolution in Omega_1

Index = 1;                                                  % Selection of the inverse solver between 1 and 2

for k_h = 2:6
for k_s = 1:size(s_list,2)
for Proposal = 1:3

% Definition of parameters

ind = ind + 1;                                              % Index of the simulation

h = 2^(-k_h);                                               % Mesh size

s = s_list(k_s);                                            % Pseudo-frequency

c = (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om);                  % Laplace transform of p(t) at s

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);                  % Boundary source


% Reconstruction of the test function

[wh, mesh] = forward_solver_method1(h, s, afun, pfun);

[Ah, submesh] = inverse_solver_1(wh, mesh, h, s, Proposal);

Ah = boundary_post_processing(Ah,mesh2);


% Computation of the exact solution on the same mesh and the error

a_exact = afun(submesh.vtx);

error_quad = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));

error_max = max(abs(Anh - A));


C(ind,:) = {Index, Proposal, h, s, Ah, submesh, error_quad, error_max};

end
end
end


% Index 2 = Resolution in Omega

Index = 2;

for k_h = 2:6
for k_s = 1:size(s_list,2)
for Extraction = 0:1

% Definition of parameters

ind = ind + 1;                                              % Index of the simulation

h = 2^(-k_h);                                               % Mesh size

s = s_list(k_s);                                            % Pseudo-frequency

c = (1-exp(-2*pi*s/om))/((1+(s/om)^2)*om);                  % Laplace transform of p(t) at s

pfun = @(X) c * (abs(X(:,2) - 1.5) < h/2);                  % Boundary source


% Solve the inverse problem through other functions

[wh, mesh] = forward_solver_method1(h, s, afun, pfun);

[Ah, submesh] = inverse_solver_2(wh, mesh, h, s, pfun, Extraction);

Ah = boundary_post_processing(Ah,mesh2);


% Calculation of the exact solution on the same mesh and the error

a_exact = afun(submesh.vtx);

error_quad = sqrt(sum((Ah-a_exact).^2))/sqrt(sum(a_exact.^2));

error_max = max(abs(Ah - A));


C(ind,:) = {Index, Extraction, h, s, Ah, submesh, error_quad, error_max};

end
end
end

% Save the data in a file

save Data_simulation_1_post.mat 'C'
