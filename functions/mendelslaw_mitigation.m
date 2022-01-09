function [R_marked,mendel_violation_with_nfp_locations] = mendelslaw_mitigation(R_marked,fp_locs)


mendel_violation_locations = mendelslaw_attack(R_marked);



idx = ~ismember(mendel_violation_locations,fp_locs,'rows');
mendel_violation_with_nfp_locations = mendel_violation_locations(idx,:);



for i = 1:size(mendel_violation_with_nfp_locations,1)
    row_id = mendel_violation_with_nfp_locations(i,1);
    col_id = ['Var' num2str(mendel_violation_with_nfp_locations(i,2)) ];
    if R_marked.(col_id)(row_id)~=1
        R_marked.(col_id)(row_id)=1;
    end
end




end