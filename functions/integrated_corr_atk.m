function R_marked_flip = integrated_corr_atk(R_marked, row_corr_pub, marginals_public,joints_public , gamma_r,gamma_l)


%% check mendel's law validation
mendel_violation_locations = mendelslaw_attack(R_marked);


%% check row-wise correlation discrepency
R_marked_content = R_marked(:,2:end).Variables;
mask  = kron( eye( size(R_marked_content ,1)/3, size(R_marked_content ,1)/3 )  ,    [0 0 1; 0 0 1; 1 1 0]);
row_corr_marked = mask.* (R_marked_content*R_marked_content');



% [row_idx,col_idx] = find( sparse(row_corr_pub) == sparse(row_corr_marked));
% 
% unchanged_rows_locs = [row_idx,col_idx];
% 
% [row_idx_zero,col_idx_zero] = find(row_corr_pub ==0);

% unfp_idx = setdiff([row_idx col_idx], [row_idx_zero col_idx_zero],   'rows','legacy'); % in A not in B 



[row_idx,col_idx] = find( sparse(row_corr_pub) ~= sparse(row_corr_marked));

fp_row = unique(  row_idx(:) );

% idx = ~ismember(unchanged_rows_locs, [row_idx_zero col_idx_zero], 'rows');% in A not in B 
% unfp_idx = unchanged_rows_locs(idx,:);
% 
% 
% 
% unfp_row = unique(unfp_idx(:,2));
% 
% fp_row = ~ismember( ([1:size(R_marked,1)])', unfp_row ,'rows'); % get all suspected rows (snp sequences)
% fp_row = find(fp_row);



%% pre-column-wise-attack (compromise entries in suspicious rows identified in the previous subsection)
R_sus_rows = R_marked(fp_row,:);

cumulative_locations = [];

tic;
[marginals_marked,joints_marked,~,~] = empirical_distributions(R_sus_rows);
toc;

att_list = R_marked.Properties.VariableNames(2:end);
for i = 1:length(att_list)
    marginal_att_i_sus = marginals_marked.(att_list{i});
    marginal_att_i_pub = marginals_public.(att_list{i});
    marginal_diff_i =    marginal_att_i_sus - marginal_att_i_pub;
    sus_entries = find(marginal_diff_i>0)-1;
    for j = 1:length(sus_entries)
        temp_idx = find( R_sus_rows.(att_list{i}) == sus_entries(j)  );
%         R_sus_rows.(att_list{i})(temp_idx) = 0;
        cumulative_locations = [cumulative_locations ;  [ fp_row(temp_idx) (i+1) * ones(length(temp_idx),1) ] ];
    end
end

% R_marked_flip = R_marked;
% R_marked_flip(fp_row,:) =  R_sus_rows;


cumulative_locations = union(cumulative_locations, mendel_violation_locations, 'rows');
R_marked_flip = flipping_attack(R_marked ,cumulative_locations);
%%

diff_thr_list = gamma_r * gamma_l /2* ones(1,2);

rounds = length(diff_thr_list);



for r  = 1:rounds
    
    
    diff_thr = diff_thr_list(r);
    
    %     R_sus_rows = R_marked_flip(fp_row,:);
    tic;
    [marginals_marked,joints_marked,~,~] = empirical_distributions( R_marked_flip  );
    toc;
    
    
    
    
    [select_row, select_col] = ...
        obtain_suspicious_row_col(joints_public, joints_marked,R_marked_flip,diff_thr);
    
    %     select_row_col1_col2 =[ [  select_row select_col(:,1) ]  ; [  select_row select_col(:,2) ]  ];
    
    if isempty(select_col), continue; end
    
    
    select_row_col1 = [ select_row   select_col(:,1)  ];
    
    idx = ~ismember(select_row_col1, cumulative_locations, 'rows');
    select_row_col1 = select_row_col1(idx,:);
    
    if isempty(select_row), continue; end;
    
    flip_location = [];
    
    for i = 1:length(fp_row)
        fp_row_i_col =  select_row_col1(  find(select_row_col1(:,1) == fp_row(i)   ),2 );
        
        if isempty(fp_row_i_col), continue; end
        
        [GC,GR] = groupcounts(  fp_row_i_col(:,1))    ;
        [val,idx] = sort(GC,'descend');
        
        %     attk_l = ceil(size(R_marked_content,2)*gamma_l/10) * 10;
        
        attk_l = floor(size(R_marked_content,2)*gamma_l);
        
        
        if length(idx)<attk_l
            fp_col =  sort(  GR(idx)  );
        else
            fp_col =  sort(  GR(idx(1:attk_l))  );
        end
        temp_loc =  [ fp_row(i) * ones(length(fp_col),1) fp_col];
        flip_location = [flip_location;temp_loc ];
    end
    
    
    flip_location = setdiff(flip_location, cumulative_locations, 'rows'); % without repetition is ok
    
    R_marked_flip = flipping_attack(R_marked_flip,flip_location);
    
    cumulative_locations = [cumulative_locations ; flip_location];
    
end

end