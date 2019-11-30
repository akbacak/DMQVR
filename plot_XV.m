clear all;
close all;
clc;

load('myDataset/hashCodes/filenames.mat');
load('myDataset/hashCodes/hashCodes_128.mat');
load('myDataset/hashCodes/targets.mat');

    N = 120;           % Number of samples in the myDataset
    data = hashCodes_128; % Binary features (Hash codes) N x NumberHasBits
    
    queryIndex = xlsread('qLabels_1.xls');  % Reads randomly choosen query pairs from excell file
    queryIndex = transpose( queryIndex ); 
    queryIndex1 = 18;       % First element of Query Pair
    queryIndex2 = 33;       % Second element of Query Pair

    l=1;
    

    
            
        q1 = data(queryIndex1,:);         % q1 & q2 are query pairs in the loop
        q2 = data(queryIndex2,:);
        q1_rep{l,:} = repmat(q1(l,:),N,1); % Make query matrix size to the same as data matrix size
        q2_rep{l,:} = repmat(q2(l,:),N,1);      
        xor_data_q1new{l,:} = xor(data, q1_rep{l,:}); % xor of data and query matrices
        xor_data_q2new{l,:} = xor(data, q2_rep{l,:});       
        hamming_dist1{l,:} = sum(xor_data_q1new{l,:},2); % sum up rows to get hamming distances
        hamming_dist2{l,:} = sum(xor_data_q2new{l,:},2);
        norm_hamming_dist1{l,:} =  hamming_dist1{l,:} / max( hamming_dist1{l,:}(:) ); % Normalize hamming  distances between 0&1
        norm_hamming_dist2{l,:} =  hamming_dist2{l,:} / max( hamming_dist2{l,:}(:) );
       
        
        X = zeros(2,N);
        X(1,:) = norm_hamming_dist1{l,:};
        X(2,:) = norm_hamming_dist2{l,:};
    
        X = (X)';
       
        maxFront = 3;
        
       [pf_idx] = pareto_fronts(X, maxFront);
        
       
       plot(X(:,1),X(:,2), 'k*','LineWidth',5); hold on 
       %h(1)= plot(X(:,1),X(:,2), 'k*','LineWidth',3); hold on          
       %h(2)= plot(pf_idx{1,1}(:,1), pf_idx{1,1}(:,2), '-b+', 'LineWidth',5); 
       %h(3)= plot(pf_idx{2,1}(:,1), pf_idx{2,1}(:,2), '-r*', 'LineWidth',5); 
       %h(4)= plot(pf_idx{3,1}(:,1), pf_idx{3,1}(:,2), '-go', 'LineWidth',5); 
       
       set(gca,'FontSize',34);  
        
       %legend(h([2 3 4]),{'First Pareto Front','Second Pareto Front','Third Pareto Front',  },'Location','southwest', 'FontSize', 34)
      
       %ylabel('d1 ', 'FontSize', 40)
       %xlabel('d2  ', 'FontSize', 40)
       
       
       
       
       
       
       
           
         
        