%*****************************************************************************************
% Edited by Enver Akbacak , 12/2018
% Marmara University, Istanbul
% akbacakk@gmail.com
% This file averages nDCG scores of all the fronts produced by
% "evaluate_by_multiple_query_pairs.m" script in the current directory
% 
% 
%
%*****************************************************************************************


 % evaluate_by_multiple_query_pairs;

 evaluate_by_multiple_query_triples;



for ll=1:l
    for jj=1:j
        [R_1(ll,jj) , C_1(ll,jj)] = size(n_rigth_DCG{ll,jj}); % Find size of each matrix in the n_rigth_DCG array
        [R_2(ll,jj) , C_2(ll,jj)] = size(n_left_DCG{ll,jj}); % Find size of each matrix in the n_left_DCG array
        
    end
end
max_1 = max(C_1(:)); % Find max size of n_rigth_DCG row
max_2 = max(C_2(:)); % Find max size of n_left_DCG row

for ll=1:l
    for jj=1:j
        n_rigth_DCG_zp{ll,jj}(1,:) =  [ n_rigth_DCG{ll,jj}(1,:) ,(zeros(1 ,max_1 - C_1(ll,jj) ))]; % Zero padding of all elements  in the n_rigth_DCG
        n_left_DCG_zp{ll,jj}(1,:)  =  [ n_left_DCG{ll,jj}(1,:)  ,(zeros(1 ,max_2 - C_2(ll,jj) ))]; % Zero padding of all elements  in the n_left_DCG
        
    end
end


n_rigth_DCG_zp_mean =[];
n_left_DCG_zp_mean =[];

for jj=1:j
    n_rigth_DCG_zp_sum  = 0;
    n_left_DCG_zp_sum  = 0;
    
    for ll=1:l
           
            n_rigth_DCG_zp_sum =  n_rigth_DCG_zp_sum  + n_rigth_DCG_zp{ll,jj}(1,:) ; % Sum of all rows for each righ front
            n_left_DCG_zp_sum =   n_left_DCG_zp_sum   + n_left_DCG_zp{ll,jj}(1,:) ;  % Sum of all rows for each left front
            
           
    end
    n_rigth_DCG_zp_mean(end+1,:) =  n_rigth_DCG_zp_sum/ll; % Mean value of each rigth front
    n_left_DCG_zp_mean(end+1,:) =  n_left_DCG_zp_sum/ll;   % Mean value of each left front
    
end

flip_n_left_DCG_zp_mean = fliplr(n_left_DCG_zp_mean);
n_DCG_mean_fronts = horzcat(flip_n_left_DCG_zp_mean, n_rigth_DCG_zp_mean);

n_DCG_mean = mean(n_DCG_mean_fronts ,1); % Mean value of first 10 fronts

plot(n_DCG_mean); 

ylabel('Mean Value of nDCG')
xlabel('Number of Retrived Items on the Fronts') 
