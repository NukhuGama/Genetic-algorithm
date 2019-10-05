%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;
clc;


%% PARAMETERS
populationSize = 40 ;
numberOfGenes = 10;
crossoverProbability = 0.7 ;
mutationProbability = 0.0525;
ParentSelectionionParameter = 0.5;
variableRange = 100   % Range for x1 and x2 between -100 up to 100 
numberOfGenerations = 250;
numberOfVariables = 2;   % Two  Variables x1 and x2
ParentSize = 10;
numberOfReplications = 2;



%% VARIABLES
fitness = zeros(populationSize, 1);

%% INITIATE POPULATION For Function1 and Function2
population1 = InitializePopulation(populationSize, numberOfGenes) ; 
population2 = InitializePopulation(populationSize, numberOfGenes) ;

%% RUN GENERATIONS For Function 1
for iGeneration = 1: numberOfGenerations
    
    %% FIND Minimum Fitness OF POPULATION
    %% Decode Chromossome of population  For Function1 and Function2
    decodedPopulation1 = DecodePopulation(population1, numberOfVariables, variableRange);
    decodedPopulation2 = DecodePopulation(population2, numberOfVariables, variableRange);
    
    % Get the minimum and maximum number from colomn as X1 and X2 
%     fitness1 = Fitness1(min(min(decodedPopulation1(:,:))), max(max(decodedPopulation1(:,:))));
%     fitness2 = Fitness2(min(min(decodedPopulation2(:,:))), max(max(decodedPopulation2(:,:))));
%   
    % Get Fitness Value from decodedPopulation 
    for i=1:length(decodedPopulation1)
        fitness1(i) = Fitness1(decodedPopulation1(i,1), decodedPopulation1(i,2));
        fitness2(i) = Fitness2(decodedPopulation2(i,1), decodedPopulation2(i,2));
    end;
    
    % Get Minimum Value of Function1 and Function2 (if fitness is maximum
    % then the value of function is minimum and if fitness is mminimum then
    % the value of function is minimum, because we pick fitness = 2^-h
    % where h is the value of function
    [minimumFunction1, bestIndividualIndex1] = max(fitness1);
    [minimumFunction2, bestIndividualIndex2] = max(fitness2);
    
    
    minimumFunction1 = -log2(minimumFunction1);
    minimumFunction2 = -log2(minimumFunction2);
    
    % Best Of X1 and X2 from Decode Chromossome 
    xBest1 = decodedPopulation1(bestIndividualIndex1,:);
    xBest2 = decodedPopulation2(bestIndividualIndex2,:);

    
    % Assign New Population of Function1 and Function2
    newPopulation1 = population1;
    newPopulation2 = population2;

    %% NEW GENERATION For Function 1
    for i = 1:ParentSize:populationSize
        %% PARENT SELECTION
        i1 = ParentSelection(fitness1,ParentSelectionionParameter,ParentSize);
        i2 = ParentSelection(fitness1,ParentSelectionionParameter,ParentSize);
        chromosome1 = population1(i1,:);
        chromosome2 = population1(i2,:);

        %% CROSS-OVER
        r = rand;
        if ( r < crossoverProbability)
            newChromosomePair = Cross(chromosome1, chromosome2);
            newPopulation1(i,:) = newChromosomePair(1,:);
            newPopulation1(i+1,:) = newChromosomePair(2,:);
        else
            newPopulation1(i,:) = chromosome1;
            newPopulation1(i+1,:) = chromosome2;
        end
    end
    
    %% NEW GENERATION For Function 2
    for i = 1:ParentSize:populationSize
        %% PARENT SELECTION
        i1 = ParentSelection(fitness2,ParentSelectionionParameter,ParentSize);
        i2 = ParentSelection(fitness2,ParentSelectionionParameter,ParentSize);
        chromosome1 = population2(i1,:);
        chromosome2 = population2(i2,:);

        %% CROSS-OVER
        r = rand;
        if ( r < crossoverProbability)
            newChromosomePair = Cross(chromosome1, chromosome2);
            newPopulation2(i,:) = newChromosomePair(1,:);
            newPopulation2(i+1,:) = newChromosomePair(2,:);
        else
            newPopulation2(i,:) = chromosome1;
            newPopulation2(i+1,:) = chromosome2;
        end
    end

    %% MUTATE
    newPopulation1 = Mutate(newPopulation1, mutationProbability);
    
    newPopulation2 = Mutate(newPopulation2, mutationProbability);
    
    %% PRESERVATION OF PREVIOUS BEST SOLUTION
    [bestChromosome1] = population1(bestIndividualIndex1,:);
    [bestChromosome2] = population2(bestIndividualIndex2,:);
    
    newPopulation1 = InsertBestIndividual(newPopulation1, bestChromosome1, numberOfReplications);
    newPopulation2 = InsertBestIndividual(newPopulation2, bestChromosome2, numberOfReplications);
        
    %% COPY THE NEW POPULATION ONTO CURRENT POPULATION
    population1 = newPopulation1;
    population2 = newPopulation2; 
    

end
 
% Change the size of vector to 4x10
bestChromosome1 = vec2mat(bestChromosome1,10);
bestChromosome2 = vec2mat(bestChromosome2,10);

% Print out
disp('Kromosom Terbaik Fungsi Pertama :')
disp(bestChromosome1)
fprintf('Nilai Minimum Fungsi Pertama =  %d\n',minimumFunction1);
fprintf('Solusi Terbaik X1 Fungsi Pertama =  %d\n',xBest1(1));
fprintf('Solusi Terbaik X2 Fungsi Pertama =  %d\n',xBest1(2)); 
disp('==================================================')
disp('Kromosom Terbaik Fungsi Kedua :')
disp(bestChromosome2)
fprintf('Nilai Minimum Fungsi Kedua =  %d\n',minimumFunction2);
fprintf('Solusi Terbaik X1 Fungsi Kedua =  %d\n',xBest2(1));
fprintf('Solusi Terbaik X2 Fungsi Kedua =  %d\n',xBest2(2)); 


