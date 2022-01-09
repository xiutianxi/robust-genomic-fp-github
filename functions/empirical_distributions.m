function [marginals,joints,joint_min,joint_max] = empirical_distributions(R)
%{
get the empirical pairwise correlation between attributes in a database
%}
atts = R.Properties.VariableNames; % the 1st att is primary key

[row_num, col_num] = size(R);

% get the marginal distributions
marginals = struct();
ins = [0 1 2];
for i = 2:col_num
    occ = [];
    for l = 1: length(ins), occ = [occ length(find( R.(atts{i}) == ins(l) ))]; end
    marginals.(atts{i}) = occ/row_num;
end


joint_min = 1;


joint_max = 0;

% get the joint distributions
joints = struct();
for i = 2:col_num-1
    ins_i = ins;
    for j  = i+1:col_num
        ins_j = ins;
        occ_joint  = [];
        for l_i = 1: length(ins_i)
            for l_j = 1: length(ins_j)
                occ_joint = [occ_joint  length( intersect(  find( R.(atts{i}) == ins_i(l_i) ), ...
                    find( R.(atts{j}) == ins_j(l_j) )   )  )];
            end
        end
                occ_joint =    (    reshape(occ_joint, length(ins_j),  length(ins_i) )   )';
                occ_joint = occ_joint/row_num;
                if joint_min>min(occ_joint(:))
                    joint_min = min(occ_joint(:));
                end
                if joint_max<max(occ_joint(:))
                    joint_max = max(occ_joint(:));
                end
        joints.(   [atts{i} '_with_' atts{j}]   ) = occ_joint;
    end
end