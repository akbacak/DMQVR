close all;
clear all;
clc;

load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset2/hashCodes/hashCodes_512.mat');
data = hashCodes_512;
load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset2/hashCodes/targets.mat');
targets = targets;
load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset2/hashCodes/filenames.mat');
filenames = filenames;


queryIndex1 = 10;
queryIndex2 = 17;
%queryIndex3 = 28;

q1 = data(queryIndex1,:); 
q2 = data(queryIndex2,:); 
%q3 = data(queryIndex3,:); 

N = length(filenames); 
q1new = repmat(q1,N,1);
q2new = repmat(q2,N,1);
%q3new = repmat(q3,N,1);

dist_1 = xor(data, q1new);
dist_2 = xor(data, q2new);
%dist_3 = xor(data, q2new);

hamming_dist1 = sum(dist_1,2);
hamming_dist2 = sum(dist_2,2);
%hamming_dist3 = sum(dist_3,2);

n_hamming_dist1 = mat2gray(hamming_dist1);
n_hamming_dist2 = mat2gray(hamming_dist2);
%n_hamming_dist3 = mat2gray(hamming_dist3);
 
 
X = zeros(2,N);
X(1,:) = n_hamming_dist1;
X(2,:) = n_hamming_dist2;
%X(3,:) = n_hamming_dist3;
X = (X)';

maxFront = 3;
[pf_idx] = pareto_fronts(X, maxFront);

q1_label = targets(queryIndex1,: ); % Label vector of Query 1
q2_label = targets(queryIndex2,: ); % Label vector of Query 2
%q3_label = targets(queryIndex3,: ); % Label vector of Query 2
        
b = or(q1_label , q2_label); % beta in the equation 7  
%b = or(b,q3_label);           %%%%%%%%%%%%%%%% 2 li sorguda bu satırı iptal et
absolute_b = nnz(b);         % Number of non-zero elements in the beta, nnz is a Matlab Func.
        
    
      for j  =1:maxFront
            R=0;
            C=0;
            Labels = targets(pf_idx{j,1}(:,3),:);   %%%%%%%%%%%%%%%% 2 li sorguda 3 olacak  
            [R , C] = size(Labels);                
               
                     
            for e = 1:R              
                MQUR_ALL(j,e) =  nnz( and(Labels(e,:) , b ) ) /  absolute_b;
                                           
            end
      end
      
for ll = 1:maxFront
    
    rtr_idx{ll,1} = [];
    rtr2_idx{ll,1} = [];

    
    [v(ll) , k(ll)] = size(MQUR_ALL(ll,:)); % k(ll) is the size of column vector of ll th front
    
    for ff = 1:k(ll)
        if  MQUR_ALL(ll,ff) == 1    
            
            rtr_idx{1,1}(end+1,:) = pf_idx{ll,1}( ff , 3); %%%%%%%%%%%%%%%% 2 li sorguda 3 olacak
            %rtr_idx{ll,:}(end+1,:) = pf_idx{ll,1}( ff , 3);
            rtr2_idx{ll,:}(end+1,:) = pf_idx{ll,1}( ff , :);
            
        end
         
    end
end

load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset2/features/features_512.mat')
features = features_512;


[M,C] = size(rtr_idx{1,1}(:,1));

f = features(rtr_idx{1,1}(:,1),:); 
f1 = features(queryIndex1,:); 
f2 = features(queryIndex2,:); 
%f3 = features(queryIndex3,:); 
f1_new = repmat(f1,M,1);
f2_new = repmat(f2,M,1);
%f3_new = repmat(f3,M,1);
dist_f1 = pdist2(f1 , f , 'euclid' );
dist_f2 = pdist2(f2 , f , 'euclid' );
%dist_f3 = pdist2(f3 , f , 'euclid' );



Y = zeros(2,M);
Y(1,:) = dist_f1;
Y(2,:) = dist_f2;
%Y(3,:) = dist_f3;
Y = (Y)';
Y2 = Y(:,1).^2 + Y(:,2).^2 ;
%Y2 = Y(:,1).^2 + Y(:,2).^2 + Y(:,3).^2 ;

Result = zeros(M,2);
Result(:,1) = Y2(:);
Result(:,2) = rtr_idx{1,1}(:,1);

final_rtr = unique(Result,'rows');

finat_rtr_idx = final_rtr(:,2);





    