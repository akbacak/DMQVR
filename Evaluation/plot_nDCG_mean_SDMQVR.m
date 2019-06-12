clear all;
clc;

x=linspace(0,350,1);



load('Evaluation/Dataset1/nDCG_mean_128');
plot(n_DCG_mean_128,'-k*');
hold on

load('Evaluation/Dataset1/nDCG_mean_256');
plot(n_DCG_mean_256,'-rx');
hold on


load('Evaluation/Dataset1/nDCG_mean_512');
plot(n_DCG_mean_512,'-b+'); 
hold on


load('Evaluation/Dataset1/nDCG_mean_1024');
plot(n_DCG_mean_1024,'-go'); 
hold on




title('nDCG averages of the first 5 pareto fronts.','FontSize', 24 )

ylabel('Mean Value of nDCGs', 'FontSize', 24)
xlabel('Number of Retrived Items on the Fronts. From middle to the left and to the rigth tails', 'FontSize', 24) 

legend({'128 bits','256 bits','512 bits','1024 bits',  },'Location','southwest', 'FontSize', 24)
