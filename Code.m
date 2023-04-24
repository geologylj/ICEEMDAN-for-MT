function ret=Code(lenchrom,bound)
%This function encodes variables into chromosomes for randomly initializing a population
% lenchrom   input : Chromosome length 
% bound      input : Value range of variables
% ret        output: The coding value of chromosomes
flag=0;
while flag==0
    pick=rand(1,length(lenchrom));
    ret=bound(:,1)'+(bound(:,2)-bound(:,1))'.*pick; 
    flag=test(lenchrom,bound,ret);     %Feasibility of testing chromosomes
end
        
