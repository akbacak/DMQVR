
%*****************************************************************************************
% Edited by Enver Akbacak, 11/2018
% Marmara University, Istanbul
% akbacakk@gmail.com
%
%
% This function computes  first 'num_fronts' pareto fronts with indexes of given input X. 
% The number of rows of X is equal to the number of images  and number  columns
% equal to the dimensions
%
% pf_idx refers to pareto points and indexes, indexes located in the last column.
%


function [pf_idx] = pareto_fronts(X, num_fronts)


pf=[];
%input=unique(X,'rows');
input=X;


front_level = num_fronts; 
 

for k=1:front_level  % Number of fronts 
     
    pop_input = [];
    temp_pf = [];
    membership = [];
    temp_input = [];
    index = [];
    temp_pf = [];
   
    
    for i = 1:size(input,1)           
        pop_input = repmat(input(i,:),size(input,1),1); 
        
        if( pop_input(1,1) < 1000 )   
            temp_input = input;
            temp_input(i,:) = ones(1,size(input,2))*1000;                  
            index = pop_input >= temp_input;    
            
            if(sum(index(:,:),2) < size(input,2)) 
                temp_pf(end+1,:)=input(i,:);  
            else
                continue;
            end
               
            
        end
                         
    end
    
  
    
    membership = ismember(X, temp_pf,'rows'); 
    membership_indexes{k,:} = find(membership); 
    pf{k,1}(:,:) = X(membership_indexes{k,:},:);
    
    pf_idx{k,1}(:,:) = [pf{k,1}(:,:) , membership_indexes{k,:}(:,1)]; 
    pf_idx{k,1}(:,:) = sortrows( pf_idx{k,1}(:,:),1);    
      
    input = setdiff(input,pf{k,1} ,'rows');
    
    
end



 
  


