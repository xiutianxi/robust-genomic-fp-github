function R_marked_flip = flipping_attack(R_marked, high_suspect)
%{
each row of high_suspect is a 3 field tuple, [sus_row sus_col
#of_occurrence], where sus_col is in the order of the original columns
ranged 2~14
high_suspect(:,2) is the original column attribute value
s_atts_ins: original states of attributes subject to fingerprinting, 13 files
%}


fp_att_list = R_marked.Properties.VariableNames;

fp_att_list = fp_att_list(2:end); % fingerprinted attribute list (exclude id and label)

flip_length = size(high_suspect,1);
all_states = [0 1 2];
for i = 1:flip_length
    row = high_suspect(i,1);
    sus_entry_bin  = dec2bin(   R_marked{ row, high_suspect(i,2) }   );
    
    if length(sus_entry_bin)==1, sus_entry_bin = ['0,' sus_entry_bin]; end
    
    sus_entry = str2num(sus_entry_bin);
    
    flip_value = min(2, bin2dec( num2str( xor(sus_entry,1) )) );
    

    R_marked{ row, high_suspect(i,2) }  = flip_value;
end
    
    R_marked_flip = R_marked;
    
    
    
end