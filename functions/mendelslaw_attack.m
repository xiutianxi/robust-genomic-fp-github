function mendel_violation_locations = mendelslaw_attack(R_marked)

[row_num, col_num] = size(R_marked);


possible_comb = [0 0 0;
    0 1 0;
    1 0 0;
    0 1 1;
    1 0 1;
    0 2 1;
    2 0 1;
    1 1 0;
    1 1 1;
    1 1 2;
    1 2 1;
    2 1 1;
    1 2 2;
    2 1 2;
    2 2 2];


% all_comb = combnk([ 0 0 0; 1 1 1; 2 2 2],3)

mendel_violation_locations = [];

for i = 1:3:row_num
    family = (R_marked([i,i+1,i+2],[2:end]).Variables)';
    
    idx = ~ismember(family,possible_comb,'rows');
    ia = find(idx);
%     [C,ia] = setdiff(family',possible_comb,'rows','legacy');
    if ~isempty(ia)
        mendel_violation_locations = [mendel_violation_locations   ;  [ kron(ones(length(ia),1),([i:1:i+2])' )  ...
            kron(ia+1,ones(3,1))]   ];
    end
end




end