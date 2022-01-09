
clc;clear;close all

all_pop = dlmread('population_w_o_snp.txt');

% parents = all_pop([1:1000],[1:101]);

parents = all_pop;

family = [];

L  = 156;

rng(100)

for i = 2:2:10000
    mother = parents(i-1,[2:end]);
    father = parents(i,[2:end]); % get rid of the index
    child = zeros(1,L);
    for j = 1:L
        if mother(j)+father(j)==0
            child(j) = 0;
            continue
        end
        if mother(j)+father(j)==1
            child(j) = rand>0.5;
            continue
        end
        if mother(j)==1 && father(j)==1
            child(j) = datasample([0 1  2],1,'Weights',[1/4 1/2 1/4]);
            continue
        end
        if mother(j)+father(j)==3
            child(j) = 1+(rand>0.5);
            continue
        end
        if mother(j)==0 && father(j)==2
            child(j) = 1;
            continue
        end
        if mother(j)==2 && father(j)==0
            child(j) = 1;
            continue
        end
        if mother(j)==2 && father(j)==2
            child(j) = 2;
            continue
        end
    end
    family = [family; mother; father; child];
    i
end

All_Pop = array2table( [([1:15000])' family]); % add the index in the first column
