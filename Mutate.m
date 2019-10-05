

function tempPopulation = Mutate(tempPopulation, mutationProbability)

indexes = rand(size(tempPopulation))<mutationProbability;                
tempPopulation(indexes) = tempPopulation(indexes)*-1+1;                     

% Deprecated - to be deleted in the next iteration
% nGenes= size(chromosome,2);
% mutatedChromosome = chromosome;
% for j = 1: nGenes
%     r= rand;
%     if (r < mutationProbability)
%         mutatedChromosome(j) = 1-chromosome(j);
%     end
%     
% end

