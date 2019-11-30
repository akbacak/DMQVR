
clear all;
close all;
clc;

 X = rand(200,2);
 
   maxFront = 3;
        
       [pf_idx] = pareto_fronts(X, maxFront);
        
        
       h(1)= plot(X(:,1),X(:,2), 'k*'); hold on
       h(2)= plot(pf_idx{1,1}(:,1), pf_idx{1,1}(:,2), '-b+'); 
       h(3)= plot(pf_idx{2,1}(:,1), pf_idx{2,1}(:,2), '-r*'); 
       h(4)= plot(pf_idx{3,1}(:,1), pf_idx{3,1}(:,2), '-go'); 
          
       set(gca,'FontSize',20);
 
       legend(h([2 3 4]),{'First Pareto Front','Second Pareto Front','Third Pareto Front',  },'Location','southwest', 'FontSize', 28)