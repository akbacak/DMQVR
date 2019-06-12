%*****************************************************************************************
% Edited by Enver Akbacak , 10/2018
% Marmara University, Istanbul
% akbacakk@gmail.com
% THIS FILE COMPUTES MQUR and  nDCG scores for multiple query pairs. 
% 
% Query pairs located in the qLabels.xls  file in the current directory )
%
%
%*****************************************************************************************

clear all;
close all;
clc;

load('myDataset/hashCodes/filenames.mat');
load('myDataset/hashCodes/hashCodes_1024.mat');
load('myDataset/hashCodes/targets.mat');

    N = 120;           % Number of samples in the myDataset
    data = hashCodes_1024; % Binary features (Hash codes) N x NumberHasBits
    
    queryIndex = xlsread('qLabels_1.xls');  % Reads randomly choosen query pairs from excell file
    queryIndex = transpose( queryIndex ); 
    queryIndex1 = queryIndex(1,:);        % First element of Query Pair
    queryIndex2 = queryIndex(2,:);        % Second element of Query Pair
    
    
    for l = 1:250                 % Number of Query Pairs
              
        q1 = data(queryIndex1,:);         % q1 & q2 are query pairs in the loop
        q2 = data(queryIndex2,:);
        q1_rep{l,:} = repmat(q1(l,:),N,1); % Make query matrix size to the same as data matrix size
        q2_rep{l,:} = repmat(q2(l,:),N,1);      
        xor_data_q1new{l,:} = xor(data, q1_rep{l,:}); % xor of data and query matrices
        xor_data_q2new{l,:} = xor(data, q2_rep{l,:});       
        hamming_dist1{l,:} = sum(xor_data_q1new{l,:},2); % sum up rows to get hamming distances
        hamming_dist2{l,:} = sum(xor_data_q2new{l,:},2);
        %norm_hamming_dist1{l,:} =  hamming_dist1{l,:} / max( hamming_dist1{l,:}(:) ); % Normalize hamming  distances between 0&1
        %norm_hamming_dist2{l,:} =  hamming_dist2{l,:} / max( hamming_dist2{l,:}(:) );
        %dist1{l,:} = mat2gray(dist1{l,:}); % Normalize hamming  distances between 0&1
        %dist2{l,:} = mat2gray(dist2{l,:});     
        
        X = zeros(2,N);
        X(1,:) = hamming_dist1{l,:};
        X(2,:) = hamming_dist2{l,:};
    
        X = (X)';
                
        maxFront = 5; 
        
        [pf_idx] = pareto_fronts(X, maxFront);   
        
        
        q1_label{l,:} = targets(queryIndex1(:,l), : ); % Label vector of Query 1
        q2_label{l,:} = targets(queryIndex2(:,l), : ); % Label vector of Query 2
        
        b{l,:} = or(q1_label{l,:} , q2_label{l,:});    % beta in the equation 7 
        absolute_b{l,:} = nnz(b{l,:});                 % Number of non-zero elements in the beta, nnz is a Matlab Func.
        
         
        for j  =1:maxFront
            
            Labels{l,j} = targets(pf_idx{j,1}(:,3),:);
                                     
            [R(l,j) , C(l,j)] = size(Labels{l,j});        
            
            
                      
            
            switch (mod(R(l,j) ,2) == 1)
                case 1
                     e_left(l,j) = (round(R(l,j) / 2) - 1) ;  
                     e_rigth(l,j)= (round(R(l,j) / 2)    ) ;
                case 0
                     e_left(l,j)  =  R(l,j) / 2 ;
                     e_rigth(l,j) =  R(l,j) / 2 ;
            end
            
            
           
            % ALL MQUR scores for each query pairs,  each front X Number of retrieved items               
            for e = 1:R(l,j)               
                MQUR_ALL{l,:}(j,e) =  nnz( and(Labels{l,j}(e,:) , b{l,:} ) ) /  absolute_b{l,:};
                                           
            end                   
            %MQUR score of the first retrieved items in each front
            MQUR_1(l,j) = MQUR_ALL{l,:}(j,1);
                        
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             rtr_idx{j,1} = [];

    
            [v(j) , k(j)] = size(MQUR_ALL{l,1}(j,:)); % k(ll) is the size of column vector of ll th front
    
             for ff = 1:k(j)
                 if  MQUR_ALL{l,1}(j,ff) == 1    
                     
                     rtr_idx{1,1}(end+1,:) = pf_idx{j,1}( ff , 3);
                     %rtr_idx{j,1}(end+1,:) = pf_idx{j,1}( ff , 3);
                     %rtr_idx{ll,:}(end+1,:) = pf_idx{ll,1}( ff , :);
                 end
             end
             
             [kk(j) zz(j)] = size( rtr_idx{j,1}(:) );
             kk = kk';
             count_rtr = sum(kk(:));
          
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
              
            
            
            
            
            
           
            
             
            G{l,j} = [];
            rigth_G{l,j} = [];
            left_G{l,j}  = [];
            i_G{l,j} =[];
            i_rigth_G{l,j} = [];
            i_left_G{l,j}  = [];
            
            rigth_DG{l,j} = [];
            left_DG{l,j}  = [];
            i_DG{l,j} = [];
            i_rigth_DG{l,j} = [];
            i_left_DG{l,j}  = [];
            
            rigth_DCG{l,j} = [];
            left_DCG{l,j}  = [];
            i_rigth_DCG{l,j} = [];
            i_left_DCG{l,j}  = [];
            
            
             for s = 1:R(l,j)     
               G{l,j}(:,s) = MQUR_ALL{l,:}(j,s);
               i_G{l,j}(:,s) = 1;
             end
             
             switch (mod(R(l,j) ,2) == 1)
                case 1                    
                    rigth_G{l,j}   =   G{l,j}(: , (e_rigth(l,j)   )  : R(l,j));
                    i_rigth_G{l,j} = i_G{l,j}(: , (e_rigth(l,j)   )  : R(l,j));
                    
                    left_G{l,j}    =   G{l,j}(: , 1 : e_left(l,j) ); 
                    left_G{l,j}    =   fliplr(left_G{l,j});
                    i_left_G{l,j}   = i_G{l,j}(: , 1 : e_left(l,j) );
             
                    rigth_DG{l,j}(:,1)   =   rigth_G{l,j}(:,1);
                    i_rigth_DG{l,j}(:,1) = i_rigth_G{l,j}(:,1);
                    
             
             
                    %left_DG{l,j}(:,1)   =   left_G{l,j}(:,1);
                    left_DG   =   left_G;
                    i_left_DG{l,j}(:,1) = i_left_G{l,j}(:,1);
             
                    for d = 2 : e_rigth(l,j)
                         rigth_DG{l,j}(:,d)   =   rigth_G{l,j}(:,d) / log2(d); %%%%%%%%%%%%%%%
                         i_rigth_DG{l,j}(:,d) = i_rigth_G{l,j}(:,d) / log2(d); 
                    end                    
                    for  f = 2 : e_left(l,j)   
                         left_DG{l,j}(:,f)   =   left_G{l,j}(:,f) / log2(f);
                         i_left_DG{l,j}(:,f) = i_left_G{l,j}(:,f) / log2(f);
                    end
                    
                    
                 case 0
                    rigth_G{l,j}   =   G{l,j}(: , (e_rigth(l,j) +1)  : R(l,j));
                    i_rigth_G{l,j} = i_G{l,j}(: , (e_rigth(l,j) +1)  : R(l,j));
                    
                    left_G{l,j}    =   G{l,j}(: , 1 : e_left(l,j) ); 
                    left_G{l,j}    =   fliplr(left_G{l,j});
                    i_left_G{l,j}   = i_G{l,j}(: , 1 : e_left(l,j) );
             
                    rigth_DG{l,j}(:,1)   =   rigth_G{l,j}(:,1);
                    i_rigth_DG{l,j}(:,1) = i_rigth_G{l,j}(:,1);
                    %rigth_DG{l,j}(:,e_rigth(l,j)) =  MQUR_ALL{l,:}(j,e_rigth(l,j));
             
             
                    left_DG{l,j}(:,1)   =   left_G{l,j}(:,1);
                    i_left_DG{l,j}(:,1) = i_left_G{l,j}(:,1);
             
             
                    for h = 2 : e_rigth(l,j)
                         rigth_DG{l,j}(:,h)   =   rigth_G{l,j}(:,h) / log2(h);
                         i_rigth_DG{l,j}(:,h) = i_rigth_G{l,j}(:,h) / log2(h); 
                    end
                    for m = 2 : e_left(l,j)
                         left_DG{l,j}(:,m)   =   left_G{l,j}(:,m) / log2(m);
                         i_left_DG{l,j}(:,m) = i_left_G{l,j}(:,m) / log2(m);
                    end
             end
             
                                    
             rigth_DCG{l,j} = cumsum(rigth_DG{l,j}); 
             i_rigth_DCG{l,j} = cumsum( i_rigth_DG{l,j});
             n_rigth_DCG{l,j} = rigth_DCG{l,j} ./ i_rigth_DCG{l,j};
                                  
           
             left_DCG{l,j} = cumsum(left_DG{l,j}); 
             i_left_DCG{l,j} = cumsum( i_left_DG{l,j});
             n_left_DCG{l,j} = left_DCG{l,j} ./ i_left_DCG{l,j};
             
             flip_n_left_DCG{l,j} = fliplr(n_left_DCG{l,j});
             n_DCG{l,j} = horzcat(flip_n_left_DCG{l,j}, n_rigth_DCG{l,j});
                 
            
           
             
        end
    end
    
 
   
    
    
    




    
      
         
