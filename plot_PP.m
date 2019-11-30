clear all;
close all;
clc;

 X = rand(200,2);
 
   maxFront = 3;
        
       [pf_idx] = pareto_fronts(X, maxFront);
        
        
 plot(X(:,1),X(:,2), 'k*'); hold on
 set(gca,'FontSize',40);
       
 legend({'First Pareto Points'},'Location','southwest', 'FontSize', 28)      
       
      