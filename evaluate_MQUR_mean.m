
%*****************************************************************************************
% Edited by Enver Akbacak , 12/2018
% Marmara University, Istanbul
% akbacakk@gmail.com
% This file averages MQUR scores of all the fronts produced by
% "evaluate_by_multiple_query_pairs.m" script in the current directory
% 
% 
%
%*****************************************************************************************


evaluate_by_multiple_query_pairs;


for ll=1:l
    
        [R(ll) , C(ll)] = size(MQUR_ALL{ll,1}); % Find size of each matrix in the MQUR_ALL array
end

max = max(C(:));  % Find max size MQUR_ALL

for ll=1:l
   for jj=1:j 
       
       % We overlapp the centers of the fronts in terms of MQUR
       switch (mod(max - C(ll,1)  ,2) == 1)
           case 1
               
                zz = (max - 1 - C(ll,1))/2;
                MQUR_ALL_zp{ll,1}(jj,:) = [ zeros(1 , zz) , MQUR_ALL{ll,1}(jj,:) , zeros(1 , zz + 1) ]; % zero pedding to the left and the rigth
              
               
           case 0
                zz = (max - C(ll,1))/2; 
                MQUR_ALL_zp{ll,1}(jj,:) = [ zeros(1 ,zz ) , MQUR_ALL{ll,1}(jj,:), zeros(1 , zz)];
              
       end
             
       
   end       
  
end


MQUR_ALL_fronts_mean =[];

for jj=1:j
    
    MQUR_ALL_zp_sum  = 0;
    
    for ll=1:l
           
           MQUR_ALL_zp_sum = MQUR_ALL_zp_sum  + MQUR_ALL_zp{ll,1}(jj,:) ; % Sum of all rows for each front
            
           
    end
    MQUR_ALL_fronts_mean(end+1,:) = MQUR_ALL_zp_sum/ll; % Mean value of each  front
  
    
end


MQUR_ALL_mean = mean(MQUR_ALL_fronts_mean ,1); % Mean value of first 10 fronts

plot( MQUR_ALL_mean); 

ylabel('Mean Value of MQUR')
xlabel('Number of Retrived Items on the Fronts') 

    
 

