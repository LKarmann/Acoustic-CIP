%% Solving Acoustic Coefficient Inverse Problem
% Lorentz Karmann, 2024

% Analysis of the results.

%% Preamble

clear;
close all;
format long;

addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openMsh;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openDom;
addpath C:\Users\loren\Dropbox\Stage\gypsilab-master\gypsilab-master\openFem;
addpath C:\Users\loren\Dropbox\Stage\Algorithmes\Functions;

load Data_simulation_1_2.mat 'C'

C = C(2:end,:);

error_mat = cell2mat(C(:,[1,3,4,5,8,9,10]));



% %% Comparison between II.1 and II.2 internal
% 
% sz = nnz(error_mat(:,1)==1);
% 
% l1 = ones([1, sz]);
% l2 = ones([1, sz]);
% l3 = ones([1, sz]);
% 
% for k = 1:sz
% 
%     ind = error_mat(:,1) == 2 & error_mat(:,2) == 1 & error_mat(:,3) == error_mat(k,3) & error_mat(:,4) == error_mat(k,4);
% 
%     [~, i1] = min([error_mat(k,5), error_mat(ind,5)]);
%     [~, i2] = min([error_mat(k,6), error_mat(ind,6)]);
%     [~, i3] = min([error_mat(k,7), error_mat(ind,7)]);
% 
%     l1(k) = i1;
%     l2(k) = i2;
%     l3(k) = i3;
% end
% 
% figure
% scatter(1:sz,l1,"filled","blue")
% hold on
% scatter(1:sz,l2,"red")
% scatter(1:sz,l3,15,"green")
% axis([0, sz+1, 0.5, 2.5]);
% legend("Quadratic", "Maximal", "Threshold");
% xlabel("Index of the simulation");
% ylabel("Index of the minimal error");
% title("Comparison between Idea II.1 and II.2 internal");
% hold off
% 
% figure
% scatter(1:sz, 100*error_mat(error_mat(:,1)==1,5), "filled")
% hold on
% scatter(1:sz, 100*repelem(error_mat(error_mat(:,1)==2 & error_mat(:,2)==1,5),3), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Relative quadratic error (in %)");
% title("Relative quadratic error between Idea II.1 and II.2 internal");
% hold off
% 
% figure
% scatter(1:sz, error_mat(error_mat(:,1)==1,7), "filled")
% hold on
% scatter(1:sz, repelem(error_mat(error_mat(:,1)==2 & error_mat(:,2)==1,7),3), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Maximal error");
% title("Maximal error after threshold cut between Idea II.1 and II.2 internal");
% hold off


% %% Comparison between II.1 and II.2 external
% 
% sz = nnz(error_mat(:,1)==1);
% 
% l1 = ones([1, sz]);
% l2 = ones([1, sz]);
% l3 = ones([1, sz]);
% 
% for k = 1:sz
% 
%     ind = error_mat(:,1) == 2 & error_mat(:,2) == 0 & error_mat(:,3) == error_mat(k,3) & error_mat(:,4) == error_mat(k,4);
% 
%     [~, i1] = min([error_mat(k,5), error_mat(ind,5)]);
%     [~, i2] = min([error_mat(k,6), error_mat(ind,6)]);
%     [~, i3] = min([error_mat(k,7), error_mat(ind,7)]);
% 
%     l1(k) = i1;
%     l2(k) = i2;
%     l3(k) = i3;
% end
% 
% figure
% scatter(1:sz,l1,"filled","blue")
% hold on
% scatter(1:sz,l2,"red")
% scatter(1:sz,l3,15,"green")
% axis([0, sz+1, 0.5, 2.5]);
% legend("Quadratic", "Maximal", "Threshold");
% xlabel("Index of the simulation");
% ylabel("Index of the minimal error");
% title("Comparison between Idea II.1 and II.2 external");
% hold off
% 
% figure
% scatter(1:sz, 100*error_mat(error_mat(:,1)==1,5), "filled")
% hold on
% scatter(1:sz, 100*repelem(error_mat(error_mat(:,1)==2 & error_mat(:,2)==0,5),3), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Relative quadratic error (in %)");
% title("Relative quadratic error between Idea II.1 and II.2 external");
% hold off
% 
% figure
% scatter(1:sz, error_mat(error_mat(:,1)==1,7), "filled")
% hold on
% scatter(1:sz, repelem(error_mat(error_mat(:,1)==2 & error_mat(:,2)==0,7),3), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Maximal error");
% title("Maximal error after threshold cut between Idea II.1 and II.2 external");
% hold off


% %% Comparison between Proposals for II.1
% 
% l = error_mat(error_mat(:,1)==1 & error_mat(:,2)==1,:);
% 
% sz = size(l,1);
% 
% l1 = ones([1, sz]);
% l2 = ones([1, sz]);
% l3 = ones([1, sz]);
% 
% for k = 1:sz
% 
%     ind2 = error_mat(:,1) == 1 & error_mat(:,2) == 2 & error_mat(:,3) == l(k,3) & error_mat(:,4) == l(k,4);
%     ind3 = error_mat(:,1) == 1 & error_mat(:,2) == 3 & error_mat(:,3) == l(k,3) & error_mat(:,4) == l(k,4);
% 
%     [~, i1] = min([l(k,5), error_mat(ind2,5), error_mat(ind3,5)]);
%     [~, i2] = min([l(k,6), error_mat(ind2,6), error_mat(ind3,6)]);
%     [~, i3] = min([l(k,7), error_mat(ind2,7), error_mat(ind3,7)]);
% 
%     l1(k) = i1;
%     l2(k) = i2;
%     l3(k) = i3;
% end
% 
% 
% figure
% scatter(1:sz,l1,"filled","blue")
% hold on
% scatter(1:sz,l2,"red")
% scatter(1:sz,l3,15,"green")
% axis([0, sz+1, 0.5, 3.5]);
% legend("Quadratic", "Maximal", "Threshold");
% xlabel("Index of the simulation");
% ylabel("Index of the minimal error");
% title("Comparison between the 3 Proposals");
% hold off


% %% Comparison between II.1 and II.2 internal, Proposal 3
% 
% l = error_mat(error_mat(:,1)==1 & error_mat(:,2)==3,:);
% 
% sz = size(l,1);
% 
% l1 = ones([1, sz]);
% l2 = ones([1, sz]);
% l3 = ones([1, sz]);
% 
% for k = 1:sz
% 
%     ind = error_mat(:,1) == 2 & error_mat(:,2) == 1 & error_mat(:,3) == error_mat(k,3) & error_mat(:,4) == error_mat(k,4);
% 
%     [~, i1] = min([l(k,5), error_mat(ind,5)]);
%     [~, i2] = min([l(k,6), error_mat(ind,6)]);
%     [~, i3] = min([l(k,7), error_mat(ind,7)]);
% 
%     l1(k) = i1;
%     l2(k) = i2;
%     l3(k) = i3;
% end
% 
% figure
% scatter(1:sz,l1,"filled","blue")
% hold on
% scatter(1:sz,l2,"red")
% scatter(1:sz,l3,15,"green")
% axis([0, sz+1, 0.5, 2.5]);
% legend("Quadratic", "Maximal", "Threshold");
% xlabel("Index of the simulation");
% ylabel("Index of the minimal error");
% title("Comparison between Idea II.1 and II.2 internal");
% hold off
% 
% figure
% scatter(1:sz, 100*l(:,5), "filled")
% hold on
% scatter(1:sz, 100*error_mat(error_mat(:,1)==2 & error_mat(:,2)==1,5), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Relative quadratic error (in %)");
% title("Relative quadratic error between Idea II.1 and II.2 internal");
% hold off
% 
% figure
% scatter(1:sz, l(:,7), "filled")
% hold on
% scatter(1:sz, error_mat(error_mat(:,1)==2 & error_mat(:,2)==1,7), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Maximal error");
% title("Maximal error after threshold cut between Idea II.1 and II.2 internal");
% hold off


% %% Comparison between II.1 and II.2 external, Proposal 3
% 
% l = error_mat(error_mat(:,1)==1 & error_mat(:,2)==3,:);
% 
% sz = size(l,1);
% 
% l1 = ones([1, sz]);
% l2 = ones([1, sz]);
% l3 = ones([1, sz]);
% 
% for k = 1:sz
% 
%     ind = error_mat(:,1) == 2 & error_mat(:,2) == 0 & error_mat(:,3) == error_mat(k,3) & error_mat(:,4) == error_mat(k,4);
% 
%     [~, i1] = min([l(k,5), error_mat(ind,5)]);
%     [~, i2] = min([l(k,6), error_mat(ind,6)]);
%     [~, i3] = min([l(k,7), error_mat(ind,7)]);
% 
%     l1(k) = i1;
%     l2(k) = i2;
%     l3(k) = i3;
% end
% 
% figure
% scatter(1:sz,l1,"filled","blue")
% hold on
% scatter(1:sz,l2,"red")
% scatter(1:sz,l3,15,"green")
% axis([0, sz+1, 0.5, 2.5]);
% legend("Quadratic", "Maximal", "Threshold");
% xlabel("Index of the simulation");
% ylabel("Index of the minimal error");
% title("Comparison between Idea II.1 and II.2 external");
% hold off
% 
% figure
% scatter(1:sz, 100*l(:,5), "filled")
% hold on
% scatter(1:sz, 100*error_mat(error_mat(:,1)==2 & error_mat(:,2)==0,5), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Relative quadratic error (in %)");
% title("Relative quadratic error between Idea II.1 and II.2 external");
% hold off
% 
% figure
% scatter(1:sz, l(:,7), "filled")
% hold on
% scatter(1:sz, error_mat(error_mat(:,1)==2 & error_mat(:,2)==0,7), "filled")
% legend("Index 1", "Index 2");
% xlabel("Index of the simulation");
% ylabel("Maximal error");
% title("Maximal error after threshold cut between Idea II.1 and II.2 external");
% hold off


% %% II.1 Proposal 3 study
% 
% l = error_mat(error_mat(:,1)==1 & error_mat(:,2)==3, 3:end);
% 
% h_list = sort(unique(l(:,1)),"descend");
% 
% figure
% hold on;
% for h = h_list'
%     scatter(l(l(:,1)==h,2), min(l(l(:,1)==h,3),2),"filled","DisplayName",sprintf("%1.6f",h));
% end
% legend('Location','northwest')
% legend("show");
% title("Relative quadratic error for different values of h and s (cut-off at 2)");
% ylabel("Relative quadratic error");
% xlabel("s");
% hold off;
% 
% 
% figure
% hold on;
% for h = h_list'
%     scatter(l(l(:,1)==h,2), l(l(:,1)==h,4),"filled","DisplayName",sprintf("%1.6f",h));
% end
% legend('Location','northwest')
% legend("show");
% title("Maximal error for different values of h and s");
% ylabel("Maximal error");
% xlabel("s");
% hold off;
% 
% 
% figure
% hold on;
% for h = h_list'
%     scatter(l(l(:,1)==h,2), l(l(:,1)==h,5),"filled","DisplayName",sprintf("%1.6f",h));
% end
% legend('Location','northwest')
% legend("show");
% title("Maximal error with a threshold for different values of h and s");
% ylabel("Maximal error with a threshold");
% xlabel("s");
% hold off;


% %% Vizualisation for various s
% 
% h = 2^(-4);
% 
% s_list = unique(error_mat(:,4));
% 
% figure
% hold on
% axis([0,1,0,1,0,4]);
% pause(5)
% 
% for s = s_list'
%     Anh = C{error_mat(:,1)==1 & error_mat(:,2)==3 & error_mat(:,3)==h & error_mat(:,4)==s,6};
%     submesh = C{error_mat(:,1)==1 & error_mat(:,2)==3 & error_mat(:,3)==h & error_mat(:,4)==s,7};
% 
%     Vh2 = fem(submesh, 'P1');
% 
%     graph(Vh2, Anh)
%     title("s = "+sprintf("%1.1f",s))
%     drawnow
% 
%     pause(3);
% end
% 
% hold off














