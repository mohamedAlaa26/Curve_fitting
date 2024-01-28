// Omar Waleed Zenhom               ID: 20206130
// Mohamed Alaa El-Deen Mohammed    ID: 20206068

#include <bits/stdc++.h>
using namespace std;

double crossoverProbability = 0.8;
double mutationProbability = 0.1;

struct chromosome {
    vector<double> genes;
    double fitness;
};


void fitnessFunction(chromosome &chromosome, int points, double x[], double y[]) {
    // Evaluate the fitness of the chromosome using mean square error
    double error = 0.0;
    for (int i = 0; i < points; i++) {
        double yCalc = 0.0;
        for (int j = 0; j < chromosome.genes.size(); j++) {
            yCalc += chromosome.genes[j] * pow(x[i], j);
        }
        error += pow(yCalc - y[i], 2.0);
    }
    chromosome.fitness = error / static_cast<double>(points);
}


void tournamentSelection(vector<chromosome> &population, int tournamentSize) {
    // Create a temporary vector to store the winners
    vector<chromosome> winners;

    // Perform tournament selection
    for (int i = 0; i < population.size(); i++) {
        chromosome winner;
        winner.fitness = numeric_limits<double>::max();  // Initialize to a large value

        // Choose n individuals randomly
        for (int j = 0; j < tournamentSize; j++) {
            int randomIndex = rand() % population.size();
            const chromosome &participant = population[randomIndex];

            // Pick the one with the highest fitness
            if (participant.fitness < winner.fitness) {
                winner = participant;
            }
        }

        // Place n copies of the winner in the mating pool
        winners.push_back(winner);
    }

    // Replace the original population with the winners
    population = winners;
}


void twoPointCrossover(chromosome &parent1, chromosome &parent2) {
    // Perform two-point crossover with the specified probability
    if (rand() / static_cast<double>(RAND_MAX) < crossoverProbability) {
        int size = parent1.genes.size();
        int point1 = rand() % size;
        int point2 = rand() % size;
        if (point1 > point2) swap(point1, point2);
        for (int i = point1; i <= point2; i++) {
            swap(parent1.genes[i], parent2.genes[i]);
        }
    }
}

void nonUniformMutation(chromosome &individual, int generation, int maxGenerations, double b) {
    for (int i = 0; i < individual.genes.size(); i++) {
        // Apply mutation with the specified probability
        if (static_cast<double>(rand()) / RAND_MAX < mutationProbability) {
            // Calculate ∆Lxi and ∆Uxi
            double deltaLxi = individual.genes[i] - (-10.0);  // LBi = -10.0
            double deltaUxi = 10.0 - individual.genes[i];    // UBi = 10.0

            // Generate r1 ∈ [0,1]
            double r1 = static_cast<double>(rand()) / RAND_MAX;

            // Determine y based on r1
            // Determine mutation direction based on r1
            double y;
            if (r1 <= 0.5) {
                y = deltaLxi;
            } else {
                y = deltaUxi;
            }

            // Calculate ∆(t,y)
            double r = static_cast<double>(rand()) / RAND_MAX;
            double tDivT = static_cast<double>(generation) / maxGenerations;
            double delta = y * (1.0 - pow(r, pow(1.0 - tDivT, b)));

            // Update the gene based on y
            if (y == deltaLxi) {
                individual.genes[i] -= delta;
            } else {
                individual.genes[i] += delta;
            }
        }
    }
}


void elitistReplacement(vector<chromosome> &population, double elitismRatio) {
    // Calculate the number of chromosomes to be copied without crossover & mutation
    int elitistCount = static_cast<int>(elitismRatio * population.size());

    // Create a temporary vector to store the best individuals
    vector<chromosome> elitistIndividuals;

    // Sort the population based on fitness
    sort(population.begin(), population.end(), [](const chromosome &a, const chromosome &b) {
        return a.fitness < b.fitness;
    });

    // Copy the best individuals to the temporary vector
    for (int i = 0; i < elitistCount; i++) {
        elitistIndividuals.push_back(population[population.size() - elitistCount + i]);
    }

    // Replace the worst individuals with the best ones from the previous generation
    for (int i = 0; i < elitistCount; i++) {
        population[i] = elitistIndividuals[i];
    }
}



int main() {
//    srand(time(NULL));
    srand(static_cast<unsigned>(time(nullptr)));

    // Open input and output files
    ifstream inputFile("curve_fitting_input.txt");
    ofstream outputFile("output.txt");

    if (!inputFile.is_open() || !outputFile.is_open()) {
        cerr << "Error opening files." << endl;
        return 1;
    }

    int datasets;
    inputFile >> datasets;

    for (int datasetIndex = 1; datasetIndex <= datasets; datasetIndex++) {
        int points, degree;
        inputFile >> points >> degree;

        double x[points], y[points];
        for (int i = 0; i < points; i++) {
            inputFile >> x[i] >> y[i];
        }

        int populationSize = 100;
        int maxGenerations = 100;

        // Initialize the best individual
        chromosome bestIndividual;
        bestIndividual.fitness = numeric_limits<double>::max(); // Initialize to a large value

        // Initialize the population with random coefficients in the range [-10, 10]
        vector<chromosome> population(populationSize);
        for (int i = 0; i < populationSize; i++) {
            population[i].genes.resize(degree + 1);
            for (int j = 0; j <= degree; j++) {
                population[i].genes[j] = (rand() / static_cast<double>(RAND_MAX)) * 20.0 - 10.0;
            }
        }

        double b = 2.0;  // You can adjust this value based on your needs

        for (int generation = 0; generation < maxGenerations; generation++) {
            // Evaluate fitness for each individual in the population
            for (int i = 0; i < populationSize; i++) {
                fitnessFunction(population[i], points, x, y);
            }

            // Update the best individual for this dataset
            for (int i = 0; i < populationSize; i++) {
                if (population[i].fitness < bestIndividual.fitness) {
                    bestIndividual = population[i];
                }
            }

            // Perform tournament selection, crossover, and mutation
            for (int i = 0; i < populationSize; i += 2) {
                tournamentSelection(population, populationSize / 2);  // Set tournament size to half of the population size
                twoPointCrossover(population[i], population[i + 1]);
                nonUniformMutation(population[i], generation, maxGenerations, b);
                nonUniformMutation(population[i + 1], generation, maxGenerations, b);
            }

            // Perform elitist replacement
            elitistReplacement(population, 0.2);

        }

        // Write the best individual for this dataset to the output file
        outputFile << "Dataset Index: " << datasetIndex << endl;
        outputFile << "Best Coefficients: ";
        for (int j = 0; j <= degree; j++) {
            outputFile << bestIndividual.genes[j] << " ";
        }
        outputFile << endl;
        outputFile << "Best Mean Square Error: " << bestIndividual.fitness << endl;
        outputFile << "---------------------------------" << endl;
    }

    // Close files
    inputFile.close();
    outputFile.close();

    return 0;
}
