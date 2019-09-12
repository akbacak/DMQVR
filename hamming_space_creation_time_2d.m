
%*****************************************************************************************
% Edited by Enver Akbacak , 9/2019
% Marmara University, Istanbul
% akbacakk@gmail.com
% This script creates pareto space 100 times for the given query pairs that
% are located in the qLabels_2.xls or qLabels_1.xls files
%  
%
%
%*****************************************************************************************





clear all;
close all;
clc;

tic

for m = 1:100

load('myDataset2/hashCodes/filenames.mat');
load('myDataset2/hashCodes/hashCodes_1024.mat');
load('myDataset2/hashCodes/targets.mat');

    N = 90;           % Number of samples in the myDataset
    data = hashCodes_1024; % Binary features (Hash codes) N x NumberHasBits
    
    queryIndex = xlsread('qLabels_2.xls');  % Reads randomly choosen query pairs from excell file
    queryIndex = transpose( queryIndex ); 
    queryIndex1 = queryIndex(1,:);        % First element of Query Pair
    queryIndex2 = queryIndex(2,:);        % Second element of Query Pair
    
    
    for l = 1:150           % Number of Query Pairs
              
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
    end
end
    
toc/100