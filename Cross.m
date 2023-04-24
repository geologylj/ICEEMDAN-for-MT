function ret=Cross(pcross,lenchrom,chrom,sizepop,bound)
%This function completes the crossover operation
% pcorss                input  : crossover probability
% lenchrom              input  : Chromosome length
% chrom     input  : chromosome complex
% sizepop               input  : chromosome complex
% ret                   output : Chromosomes after crossing
 for i=1:sizepop  
% In each round of for loop, there may be a crossover operation, 
% where chromosomes are randomly selected and the crossover position is also randomly selected. 
% However, whether to perform crossover operation in this round of for loop is determined by the crossover probability (continue control)
     % Randomly select two chromosomes for crossing
     pick=rand(1,2);
     while prod(pick)==0
         pick=rand(1,2);
     end
     index=ceil(pick.*sizepop);
     % The probability of crossing determines whether to cross or not
     pick=rand;
     while pick==0
         pick=rand;
     end
     if pick>pcross
         continue;
     end
     flag=0;
     while flag==0
         pick=rand;
         while pick==0
             pick=rand;
         end
         pos=ceil(pick.*sum(lenchrom)); 
         pick=rand; 
         v1=chrom(index(1),pos);
         v2=chrom(index(2),pos);
         chrom(index(1),pos)=pick*v2+(1-pick)*v1;
         chrom(index(2),pos)=pick*v1+(1-pick)*v2; 
         flag1=test(lenchrom,bound,chrom(index(1),:)); 
         flag2=test(lenchrom,bound,chrom(index(2),:)); 
         if   flag1*flag2==0
             flag=0;
         else flag=1;
         end   
     end
 end
ret=chrom;