%% Function created by Alexis Perez Bellido 5/2014. The purpose of this function is to create the different
%  levels of an experiment and combine them. It also creates a vector with
%  a random trial order.

function [d,e] = combinator(x, NRepetitions) % [RandTrials allCombinations] = combinator([2,3,4],5)   %para llamarla

factors = x;  % you only have to change this line, the rest of the code prepares the combinations automaticallty

nfactors= size(factors);
allCombinations = zeros( 1,nfactors(2)); % Here I create the initial array

allCombinations = repmat(allCombinations ,factors(1),1);
allCombinations(:,1) = [1:factors(1)];

for h =  2 : nfactors(2)

allCombinations = repmat(allCombinations ,factors(h),1);
a=1;
len = length(allCombinations);
cocient =  len/factors(h);

for i = 1 : len
     allCombinations(i,h)=a;  
    if rem(i,cocient)==0
    a=a +1; 
    end; 
    
end;
end;

sizeComb= length(allCombinations);
NTrials=NRepetitions*sizeComb;
allCombinations=repmat(allCombinations,NRepetitions,1);
%generic variable to randomize inside the allCombination
sizeComb2=size(allCombinations);
f=1:sizeComb2(1);
ff=f(randperm(length(f)));  % This will serve to select the corresponding arrow during the experiment.
d=ff;
e=allCombinations;
 end



