

function [R_marked,fp,fp_locs] = vanilla_insert_fingerprint(R, gamma_r,gamma_l,secretKey,sp_id)
%{
obtain the marked database for a service provider (sp_id), secretKey is the
        secret key of the database owner
epsilon: mark one of the last significant epsilon bit
gamma: mark every gamma tupple
%}

[row_num,col_num] = size(R);
fp =  sp_id_fingerprint_generate(secretKey, sp_id);
L = length(fp);
Start = 1;
Stop = row_num*col_num;

R_marked = R;

fp_locs = [];

cnt = 0;




for row  = 1:row_num
    
    
    primary_key_att_value = R{row,1};
    seed_row = [double(secretKey) primary_key_att_value];
    rng(sum(seed_row));
    %     rnd1 = datasample([1:Stop],1);
    rnd1 = floor(rand*Stop);
    if ~mod(rnd1,  floor(1/gamma_r) ) % if Modulo is 0, then fingerprint this SNP sequence


        for p = 2:col_num  % the first attribute is the primary key, which is not fped
            
            
            seed_col = [double(secretKey) primary_key_att_value p];
            
            rng(sum(seed_col));
            %             rnd_seq = datasample([1:Stop],4,'Replace',false);
            rnd_seq = floor(rand(1,4)*Stop);
            if ~mod(rnd_seq(1),  floor(1/gamma_l) ) % if Modulo is 0, then fingerprint this SNP
%                 display(p)
                fp_locs  = [fp_locs ; [row p]];
                mask_bit =  mod(rnd_seq(2),2);
                fp_index =  mod(rnd_seq(3),L)+1;
                mark_bit = xor(mask_bit,fp(fp_index));
                t = mod(rnd_seq(4),2)+1;
                att_value_bin  = dec2bin( R_marked{primary_key_att_value,p});
                if length(att_value_bin)==1, att_value_bin = ['0' att_value_bin]; end
                att_value_bin(3-t) = int2str(mark_bit);
                att_value_update = bin2dec(att_value_bin);
                R_marked{primary_key_att_value,p} = min(2,att_value_update); % post-processing to eliminate snp = 3
%                 R_marked{primary_key_att_value,p} =  att_value_update;
            end
        end
    end

end


end



