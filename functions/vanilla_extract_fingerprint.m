function [f_detect,f_vote0,f_vote1] = vanilla_extract_fingerprint(R_marked, gamma_r,gamma_l,secretKey)
%{
obtain the marked database for a service provider (sp_id), secretKey is the
        secret key of the database owner
epsilon: mark one of the last significant epsilon bit
gamma: mark every gamma tupple
%}

L = 128;
f_vote0 = zeros(1,L) ;
f_vote1 = zeros(1,L) ;

[row_num,col_num] = size(R_marked);


Start = 1;
Stop = row_num*col_num;




for row  = 1:row_num
    primary_key_att_value = R_marked{row,1};
    seed_row = [double(secretKey) primary_key_att_value];
    rng(sum(seed_row));
%     rnd1 = datasample([1:Stop],1);
    rnd1 = floor(rand*Stop);
    if ~mod(rnd1,  floor(1/gamma_r) ) % if Modulo is 0, then fingerprint this SNP sequence
%         display(row)
        
        for p = 2:col_num  % the first attribute is the primary key, which is not fped
            
            
            seed_col = [double(secretKey) primary_key_att_value p];
            
            rng(sum(seed_col));
%             rnd_seq = datasample([1:Stop],4,'Replace',false);
            rnd_seq = floor(rand(1,4)*Stop);
            if ~mod(rnd_seq(1),  floor(1/gamma_l) ) % if Modulo is 0, then fingerprint this SNP
                mask_bit =  mod(rnd_seq(2),2);
                fp_index =  mod(rnd_seq(3),L)+1;
                %                 mark_bit = xor(mask_bit,fp(fp_index));
                t = mod(rnd_seq(4),2)+1;
                att_value_bin  = dec2bin( R_marked{primary_key_att_value,p});
                if length(att_value_bin)==1, att_value_bin = ['0' att_value_bin]; end
                mark_bit = bin2dec( att_value_bin(3-t) );
                f_l = xor(mark_bit,mask_bit);
                if f_l==0
                    f_vote0(fp_index) = f_vote0(fp_index) +1;
                else
                    f_vote1(fp_index) = f_vote1(fp_index) +1;
                end
            end
        end
    end
    
end




f_detect =  double(  (  f_vote1./(f_vote0+f_vote1)  )>0.5   );


f_detect(  find(  (f_vote1==f_vote0)  ==1  )   ) = nan;
end