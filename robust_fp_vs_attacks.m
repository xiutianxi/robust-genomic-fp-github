

clear;clc;close all;
addpath('saved data')
addpath('functions')

load 'All_Pop.mat'
R = All_Pop;
R = All_Pop([1:1500],:);

%% get ground truth column- and row-wise correlations

tic;
[marginals_public,joints_public,joint_min,joint_max] = empirical_distributions(R);
toc;

R_content = R(:,2:end).Variables;
mask  = kron( eye( size(R ,1)/3, size(R ,1)/3 )  ,    [0 0 1; 0 0 1; 1 1 0]);
row_corr_pub = mask.* (R_content*R_content');

%% fingerprint insertion
secretKey = 'i am the database owner';
sp_id = 100;


Gamma_R = [0.05:0.05:0.3];

Gamma_L = [0.1:0.02:0.2];

comp_fp_count = zeros(length(Gamma_R), length(Gamma_L));

chg_entry_count = zeros(length(Gamma_R), length(Gamma_L));

for g_r_idx = 1:length(Gamma_R)
    
    gamma_r  = Gamma_R(g_r_idx);
    
    for g_l_idx = 1:length(Gamma_L)
        
        gamma_l = Gamma_L(g_l_idx);
        tic;
        [R_marked,fp,fp_locs]  = vanilla_insert_fingerprint(R, gamma_r,gamma_l,secretKey,sp_id);
        toc;
        L = length(fp);
    
        
        
        R_marked = integrated_corr_mitigation(R, R_marked,fp_locs);

        %% integrated correlatin attack
        R_marked_flip = integrated_corr_atk2(R_marked, row_corr_pub, marginals_public,joints_public , gamma_r,gamma_l);
        tic;
        [f_detect,f_vote0,f_vote1] = vanilla_extract_fingerprint(R_marked_flip, gamma_r,gamma_l,secretKey);
        toc;
        comp_fp_count(g_r_idx, g_l_idx) = 1- length( intersect(find(fp==f_detect),find(isnan(f_detect)==0)) )   /  (  L-sum(isnan(f_detect))  )
        % xxx = [fp;f_detect;f_vote0;f_vote1];
        chg_entry_count(g_r_idx, g_l_idx) = sum(sum(R_marked.Variables~=R_marked_flip.Variables))/(prod(size(R_marked)))
    end
end