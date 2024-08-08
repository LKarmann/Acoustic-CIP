%% Analysis of quantitative reconstruction methods for solution to acoustic coefficient inverse problem
% Lorentz Karmann, 2024

% Analysis of the results from Simulation 1 (see Appendix).


%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions;


%% Load the data

% Without boundary post-processing

% Test 1: 1 inclusion
load C:\Users\loren\Dropbox\Stage\Algorithmes\Github\Data_simulation_1\Data_simulation_1.mat 'C'

% Test 2: 2 inclusions
% load C:\Users\loren\Dropbox\Stage\Algorithmes\Github\Data_simulation_1\Data_simulation_1_2.mat 'C'

% Extraction of the numerical values
error_mat_1 = cell2mat(C(2:end,[1,2,3,4,7,8,9]));


% With boundary post-processing

% Test 1: 1 inclusion
load C:\Users\loren\Dropbox\Stage\Algorithmes\Github\Data_simulation_1\Data_simulation_1_post.mat 'C'

% Test 2: 2 inclusions
% load C:\Users\loren\Dropbox\Stage\Algorithmes\Github\Data_simulation_1\Data_simulation_1_2_post.mat 'C'

% Extraction of the numerical values
error_mat_2 = cell2mat(C(2:end,[1,2,3,4,7,8]));



%% Comparison between Proposals for Index = 1 (without boundary post-processing)

% Extraction of the data with Index = 1 and Proposal = 1
l = error_mat_1(error_mat_1(:,1)==1 & error_mat_1(:,2)==1,:);

sz = size(l,1);

l1 = ones([1, sz]);
l2 = ones([1, sz]);

for k = 1:sz

    ind2 = error_mat_1(:,1) == 1 & error_mat_1(:,2) == 2 & error_mat_1(:,3) == l(k,3) & error_mat_1(:,4) == l(k,4);
    ind3 = error_mat_1(:,1) == 1 & error_mat_1(:,2) == 3 & error_mat_1(:,3) == l(k,3) & error_mat_1(:,4) == l(k,4);

    [~, i1] = min([l(k,5), error_mat_1(ind2,5), error_mat_1(ind3,5)]);
    [~, i2] = min([l(k,6), error_mat_1(ind2,6), error_mat_1(ind3,6)]);

    l1(k) = i1;
    l2(k) = i2;
end


figure
scatter(1:sz,l1,"filled","blue")
hold on
scatter(1:sz,l2,"red")
axis([0, sz+1, 0.5, 3.5]);
legend("Quadratic", "Maximal");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2,3]);
ylabel("Index of the proposal",'Interpreter','latex','FontSize',18);
title("Comparison of the proposals",'Interpreter','latex','FontSize',18);
hold off


%% Comparison between Proposals for Index = 1 (with boundary post-processing)

% Extraction of the data with Index = 1 and Proposal = 1
l = error_mat_2(error_mat_2(:,1)==1 & error_mat_2(:,2)==1,:);

sz = size(l,1);

l1 = ones([1, sz]);
l2 = ones([1, sz]);

for k = 1:sz

    ind2 = error_mat_2(:,1) == 1 & error_mat_2(:,2) == 2 & error_mat_2(:,3) == l(k,3) & error_mat_2(:,4) == l(k,4);
    ind3 = error_mat_2(:,1) == 1 & error_mat_2(:,2) == 3 & error_mat_2(:,3) == l(k,3) & error_mat_2(:,4) == l(k,4);

    [~, i1] = min([l(k,5), error_mat_2(ind2,5), error_mat_2(ind3,5)]);
    [~, i2] = min([l(k,6), error_mat_2(ind2,6), error_mat_2(ind3,6)]);

    l1(k) = i1;
    l2(k) = i2;
end


figure
scatter(1:sz,l1,"filled","blue")
hold on
scatter(1:sz,l2,"red")
axis([0, sz+1, 0.5, 3.5]);
legend("Quadratic", "Maximal");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2,3]);
ylabel("Index of the proposal",'Interpreter','latex','FontSize',18);
title("Comparison of the proposals",'Interpreter','latex','FontSize',18);
hold off


%% Comparison between Index 1 and Index 2 for Proposal 3 (internal, without boundary post-processing)

% Extraction of the data with Index = 1 and Proposal = 3
l = error_mat_1(error_mat_1(:,1)==1 & error_mat_1(:,2)==3,:);

sz = size(l,1);

l1 = ones([1, sz]);
l2 = ones([1, sz]);

for k = 1:sz

    ind = error_mat_1(:,1) == 2 & error_mat_1(:,2) == 1 & error_mat_1(:,3) == l(k,3) & error_mat_1(:,4) == l(k,4);

    [~, i1] = min([l(k,5), error_mat_1(ind,5)]);
    [~, i2] = min([l(k,6), error_mat_1(ind,6)]);

    l1(k) = i1;
    l2(k) = i2;
end

figure
scatter(1:sz,l1,"filled","blue")
hold on
scatter(1:sz,l2,"red")
axis([0, sz+1, 0.5, 2.5]);
legend("Quadratic", "Maximal");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2]);
ylabel("Index of the inverse solver",'Interpreter','latex','FontSize',18);
title("Comparison between the inverse methods",'Interpreter','latex','FontSize',18);
hold off


%% Comparison between Index 1 and Index 2 for Proposal 3 (internal, with boundary post-processing)

% Extraction of the data with Index = 1 and Proposal = 3
l = error_mat_2(error_mat_2(:,1)==1 & error_mat_2(:,2)==3,:);

sz = size(l,1);

l1 = ones([1, sz]);
l2 = ones([1, sz]);

for k = 1:sz

    ind = error_mat_2(:,1) == 2 & error_mat_2(:,2) == 1 & error_mat_2(:,3) == l(k,3) & error_mat_2(:,4) == l(k,4);

    [~, i1] = min([l(k,5), error_mat_2(ind,5)]);
    [~, i2] = min([l(k,6), error_mat_2(ind,6)]);

    l1(k) = i1;
    l2(k) = i2;
end

figure
scatter(1:sz,l1,"filled","blue")
hold on
scatter(1:sz,l2,"red")
axis([0, sz+1, 0.5, 2.5]);
legend("Quadratic", "Maximal");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2]);
ylabel("Index of the inverse solver",'Interpreter','latex','FontSize',18);
title("Comparison between the inverse methods",'Interpreter','latex','FontSize',18);
hold off


ll = error_mat_2(error_mat_2(:,1)==2 & error_mat_2(:,2)==1,:);

figure
scatter(1:sz,min(l(:,5),2),"filled","blue")
hold on
scatter(1:sz,min(ll(:,5),2),"filled","red")
legend("Index 1", "Index 2");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2]);
ylabel("$\frac{||a_h - a||_2}{||a||_2}$",'Interpreter','latex','FontSize',18);
title("Comparison between the inverse methods",'Interpreter','latex','FontSize',18);
hold off

%% Comparison between Index 1 and Index 2 for Proposal 3 (external, without boundary post-processing)

% Extraction of the data with Index = 1 and Proposal = 3
l = error_mat_1(error_mat_1(:,1)==1 & error_mat_1(:,2)==3,:);

sz = size(l,1);

l1 = ones([1, sz]);
l2 = ones([1, sz]);

for k = 1:sz

    ind = error_mat_1(:,1) == 2 & error_mat_1(:,2) == 0 & error_mat_1(:,3) == l(k,3) & error_mat_1(:,4) == l(k,4);

    [~, i1] = min([l(k,5), error_mat_1(ind,5)]);
    [~, i2] = min([l(k,6), error_mat_1(ind,6)]);

    l1(k) = i1;
    l2(k) = i2;
end

figure
scatter(1:sz,l1,"filled","blue")
hold on
scatter(1:sz,l2,"red")
axis([0, sz+1, 0.5, 2.5]);
legend("Quadratic", "Maximal");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2]);
ylabel("Index of the inverse solver",'Interpreter','latex','FontSize',18);
title("Comparison between the inverse methods",'Interpreter','latex','FontSize',18);
hold off


%% Comparison between Index 1 and Index 2 for Proposal 3 (external, with boundary post-processing)

% Extraction of the data with Index = 1 and Proposal = 3
l = error_mat_2(error_mat_2(:,1)==1 & error_mat_2(:,2)==3,:);

sz = size(l,1);

l1 = ones([1, sz]);
l2 = ones([1, sz]);

for k = 1:sz

    ind = error_mat_2(:,1) == 2 & error_mat_2(:,2) == 0 & error_mat_2(:,3) == l(k,3) & error_mat_2(:,4) == l(k,4);

    [~, i1] = min([l(k,5), error_mat_2(ind,5)]);
    [~, i2] = min([l(k,6), error_mat_2(ind,6)]);

    l1(k) = i1;
    l2(k) = i2;
end

figure
scatter(1:sz,l1,"filled","blue")
hold on
scatter(1:sz,l2,"red")
axis([0, sz+1, 0.5, 2.5]);
legend("Quadratic", "Maximal");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2]);
ylabel("Index of the inverse solver",'Interpreter','latex','FontSize',18);
title("Comparison between the inverse methods",'Interpreter','latex','FontSize',18);
hold off


%% Efficiency of the boundary post-processing

% Extraction of the data with Index = 1 and Proposal = 3
l = error_mat_1(error_mat_1(:,1)==1 & error_mat_1(:,2)==3,:);

sz = size(l,1);

l1 = ones([1, sz]);
l2 = ones([1, sz]);

for k = 1:sz

    ind = error_mat_2(:,1) == 1 & error_mat_2(:,2) == 3 & error_mat_2(:,3) == l(k,3) & error_mat_2(:,4) == l(k,4);

    [~, i1] = min([l(k,5), error_mat_2(ind,5)]);
    [~, i2] = min([l(k,6), error_mat_2(ind,6)]);

    l1(k) = i1;
    l2(k) = i2;
end

figure
scatter(1:sz,l1,"filled","blue")
hold on
scatter(1:sz,l2,"red")
axis([0, sz+1, 0.5, 2.5]);
legend("Quadratic", "Maximal");
xlabel("Index of the simulation",'Interpreter','latex','FontSize',18);
yticks([1,2]);
yticklabels(["Without", "With"]);
ylabel("Boundary post-processing",'Interpreter','latex','FontSize',18);
title("Efficiency of the boundary post-processing",'Interpreter','latex','FontSize',18);
hold off


%% Analysis of h and s

l = error_mat_2(error_mat_2(:,1)==1 & error_mat_2(:,2)==3, 3:end);

h_list = sort(unique(l(:,1)),"descend");

figure
hold on;
for k = 2:6
    scatter(l(l(:,1)==h_list(k-1),2), min(l(l(:,1)==h_list(k-1),3),2),"filled","DisplayName","h = 2^{-"+sprintf("%1.0f",k)+"}");
end
legend('Location','northwest')
legend("show");
title("Error analysis for various $h$ and $s$ (cut-off at 2)",'Interpreter','latex','FontSize',18);
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
xlabel("$s$",'Interpreter','latex','FontSize',18);
hold off;


figure
hold on;
for k = 2:6
    scatter(l(l(:,1)==h_list(k-1),2), l(l(:,1)==h_list(k-1),4),"filled","DisplayName","h = 2^{-"+sprintf("%1.0f",k)+"}");
end
legend('Location','northwest')
legend("show");
title("Error analysis for various $h$ and $s$",'Interpreter','latex','FontSize',18);
ylabel("$||a_h - a||_{\infty}$",'Interpreter','latex','FontSize',18);
xlabel("s",'Interpreter','latex','FontSize',18);
hold off;


s = 4.5;

x_bar = mean(log2(h_list));
y_bar = mean(log2(l(l(:,2)==s,3)));
a = sum((log2(h_list)-x_bar).*(log2(l(l(:,2)==s,3))-y_bar))/sum((log2(h_list)-x_bar).^2);
b = y_bar-a*x_bar;

figure
scatter(log2(h_list),log2(l(l(:,2)==s,3)),'filled')
hold on
plot(log2(h_list), a*log2(h_list)+b)
hold off
legend(["","$y = ax + b$ with $a = $"+sprintf('%1.4f',a)+", $b = $"+sprintf('%1.4f',b)],'Interpreter','latex')
title("Convergence in $h$ at $s =$ "+sprintf('%1.1f',s),'Interpreter','latex','FontSize',18);
ylabel("$\log_2 \frac{||a_h - a||_2}{||a||_2}$",'Interpreter','latex','FontSize',18);
xlabel("$\log_2 h$",'Interpreter','latex','FontSize',18);


s = 7;

x_bar = mean(log2(h_list));
y_bar = mean(log2(l(l(:,2)==s,4)));
a = sum((log2(h_list)-x_bar).*(log2(l(l(:,2)==s,4))-y_bar))/sum((log2(h_list)-x_bar).^2);
b = y_bar-a*x_bar;

figure
scatter(log2(h_list),log2(l(l(:,2)==s,4)),'filled')
hold on
plot(log2(h_list), a*log2(h_list)+b)
hold off
legend(["","$y = ax + b$ with $a = $"+sprintf('%1.4f',a)+", $b = $"+sprintf('%1.4f',b)],'Interpreter','latex')
title("Convergence in $h$ at $s =$ "+sprintf('%1.1f',s),'Interpreter','latex','FontSize',18);
ylabel("$\log_2 ||a_h - a||_{\infty}$",'Interpreter','latex','FontSize',18);
xlabel("$\log_2 h$",'Interpreter','latex','FontSize',18);



%% Analysis of h and s (Index = 2)

l = error_mat_2(error_mat_2(:,1)==2 & error_mat_2(:,2)==1, 3:end);

h_list = sort(unique(l(:,1)),"descend");

figure
hold on;
for k = 2:6
    scatter(l(l(:,1)==h_list(k-1),2), min(l(l(:,1)==h_list(k-1),3),2),"filled","DisplayName","h = 2^{-"+sprintf("%1.0f",k)+"}");
end
legend('Location','northwest')
legend("show");
title("Error analysis for various $h$ and $s$ (cut-off at 2)",'Interpreter','latex','FontSize',18);
ylabel("$\frac{||a_h - a||_2}{||a||_2}$","Interpreter","latex","FontSize",18);
xlabel("$s$",'Interpreter','latex','FontSize',18);
hold off;


figure
hold on;
for k = 2:6
    scatter(l(l(:,1)==h_list(k-1),2), l(l(:,1)==h_list(k-1),4),"filled","DisplayName","h = 2^{-"+sprintf("%1.0f",k)+"}");
end
legend('Location','northwest')
legend("show");
title("Error analysis for various $h$ and $s$",'Interpreter','latex','FontSize',18);
ylabel("$||a_h - a||_{\infty}$",'Interpreter','latex','FontSize',18);
xlabel("s",'Interpreter','latex','FontSize',18);
hold off;


s = 4.5;

x_bar = mean(log2(h_list));
y_bar = mean(log2(l(l(:,2)==s,3)));
a = sum((log2(h_list)-x_bar).*(log2(l(l(:,2)==s,3))-y_bar))/sum((log2(h_list)-x_bar).^2);
b = y_bar-a*x_bar;

figure
scatter(log2(h_list),log2(l(l(:,2)==s,3)),'filled')
hold on
plot(log2(h_list), a*log2(h_list)+b)
hold off
legend(["","$y = ax + b$ with $a = $"+sprintf('%1.4f',a)+", $b = $"+sprintf('%1.4f',b)],'Interpreter','latex')
title("Convergence in $h$ at $s =$ "+sprintf('%1.1f',s),'Interpreter','latex','FontSize',18);
ylabel("$\log_2 \frac{||a_h - a||_2}{||a||_2}$",'Interpreter','latex','FontSize',18);
xlabel("$\log_2 h$",'Interpreter','latex','FontSize',18);


s = 7;

x_bar = mean(log2(h_list));
y_bar = mean(log2(l(l(:,2)==s,4)));
a = sum((log2(h_list)-x_bar).*(log2(l(l(:,2)==s,4))-y_bar))/sum((log2(h_list)-x_bar).^2);
b = y_bar-a*x_bar;

figure
scatter(log2(h_list),log2(l(l(:,2)==s,4)),'filled')
hold on
plot(log2(h_list), a*log2(h_list)+b)
hold off
legend(["","$y = ax + b$ with $a = $"+sprintf('%1.4f',a)+", $b = $"+sprintf('%1.4f',b)],'Interpreter','latex')
title("Convergence in $h$ at $s =$ "+sprintf('%1.1f',s),'Interpreter','latex','FontSize',18);
ylabel("$\log_2 ||a_h - a||_{\infty}$",'Interpreter','latex','FontSize',18);
xlabel("$\log_2 h$",'Interpreter','latex','FontSize',18);