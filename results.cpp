#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;
/*
step,fitness,genome
0,0.11423977257274615,"[True, True, False, False, True, True, True, False, False, False, False, True, False, False, True, False, True, True, False, False, False, False, True, False, True, True, True, False, True, True, True, False, True, True, False, False, True, True, True, False]"
0,0.012605244161898486,"[True, False, True, False, False, True, False, True, False, True, False, False, True, True, True, False, True, True, True, True, False, True, True, True, True, False, False, False, False, True, False, True, False, True, True, False, True, True, True, False]"
*/
const int populations[] = {25, 50, 75, 100};
const float mutations[] = {0, 0.01, 0.03, 0.05};
const float probabilities[] = {0, 0.1, 0.3, 0.5};
const int t_sizes[] = {2, 3, 4, 5};

// 5120 Experiments, 30 generations per experiment
// N,p_m,p_c,tournament_size,iteration,generation,average_fitness,best_fitness,best_genome,solution_found,num_solutions_found,diversity_metric
class Experiment {
    public:
        Experiment(
            int N, 
            float p_m, 
            float p_c, 
            int tournament, 
            int iteration,
            int generation,
            long double averageFitness,
            pair<long double, string> bestPair,
            int solutions,
            vector <string> genomes
        );

        void Print();
        double calculateDivMetric(vector <string>& genomes);

        int N;
        float p_m;
        float p_c;
        int tournament;
        int iteration;
        int generation;
        long double averageFitness;
        long double bestFitness;
        string bestGenome;
        bool solutionFound;
        int solutions;
        double diversityMetric;
};

vector <Experiment> rows;
string convertToBinaryString(const string& input);
double euclideanDistance(const string& g1, const string& g2);
string formatFloat(float value);
string getFileName(int N, float p_m, float p_c, int size, int iteration);
void readCSVFile(const string& fileName, int pop, float mutation, float probability, int t_size, int iteration);
void writeCSVFile();

int main()
{
    //vector <Experiment> rows;
    // Loop through all fileNames (5120 iterations)
        // For each file:
            // create 30 Experiment instances
            // skip lines where step is a string
            // If not, increment a counter for each line
            // If line % population = 0, create an Experiment instance
    int n = 25; float m = 0.01; float p = 0.1; int size = 2;

    for (int pop: populations) {
        for (float mutation: mutations) {
            for (float probability: probabilities) {
                for (int t_size: t_sizes) {
                    for (int iteration = 1; iteration <= 20; iteration++) {
                        string file = getFileName(pop, mutation, probability, t_size, iteration);
                        readCSVFile(file, pop, mutation, probability, t_size, iteration);
                    }
                }
            }
        }

        cout << "Population: " << pop << endl;
    }

    cout << "Rows: " << rows.size() << endl;
    writeCSVFile();

    return 0;
}

string convertToBinaryString(const string& input) 
{
    std::string result;
    std::stringstream ss(input);
    std::string word;

    // Remove the first and last brackets
    std::string cleaned = input.substr(1, input.size() - 2);  

    std::stringstream stream(cleaned);
    while (std::getline(stream, word, ',')) {  
        // Trim spaces
        while (!word.empty() && (word.front() == ' ' || word.front() == '[')) 
            word.erase(word.begin());
        while (!word.empty() && (word.back() == ']' || word.back() == ' ')) 
            word.pop_back();

        if (word == "True") {
            result += '1';
        } else if (word == "False") {
            result += '0';
        }
    }

    return result;
}

double euclideanDistance(const string& g1, const string& g2) 
{
    double sum = 0.0;
    for (size_t i = 0; i < g1.size(); i++) {
        double diff = (g1[i] - g2[i]);  // Convert '0' or '1' to int and compute difference
        sum += diff * diff;  // Square the difference
    }
    return sqrt(sum);  // Take square root
}

string formatFloat(float value) 
{
    ostringstream oss;
    oss << std::defaultfloat << value; 
    return oss.str();
}

string getFileName(int N, float p_m, float p_c, int size, int iteration) 
{
    return "./experiments/" + to_string(N) + "_" + 
           formatFloat(p_m) + "_" + formatFloat(p_c) + "_" + 
           to_string(size) + "_" + to_string(iteration) + ".csv";
}

void readCSVFile(const string& fileName, int pop, float mutation, float probability, int t_size, int iteration) 
{
    string line;
    ifstream fin(fileName);
    if(!fin.is_open()){
        cout << "file " << fileName << " failed to open."<<endl;
        return;
    }

    string cell;
    int counter = 0;
    long double totalFitness = 0;
    pair <long double, string> bestGenome = {-1, ""};
    vector <string> genomes;
    int solutions = 0;

    while (getline(fin, line)) {
        if (line[4] != ',') {
            istringstream ss(line);
            int step; double fitness; string genome;
            getline(ss, cell, ',');
            step = stoi(cell);

            getline(ss, cell, ',');
            fitness = stold(cell);

            getline(ss, cell);
            genome = convertToBinaryString(cell);

            genomes.push_back(genome);
            totalFitness += fitness;
            if (fitness >= 1) solutions++;
            if (fitness > bestGenome.first) bestGenome = {fitness, genome};

            counter++;
            if (counter % pop == 0) {
                Experiment e(pop, mutation, probability, t_size, iteration, counter/pop, totalFitness/pop, bestGenome, solutions, genomes);
                rows.push_back(e);
                bestGenome = {-1, ""};
                solutions = 0;
                totalFitness = 0;
                genomes.clear();
            }
        }
    }

    fin.close();
}

void writeCSVFile()
{
    ofstream file("Final.csv");
    if (!file.is_open()) {
        cout << "Error opening file!" << std::endl;
        return;
    }

    for (Experiment e: rows) {
        std::ostringstream avgFit, bestFit;
        avgFit << std::fixed << std::setprecision(15) << e.averageFitness;
        bestFit << std::fixed << std::setprecision(15) << e.bestFitness;

        file << e.N << "," << e.p_m << "," << e.p_c << "," << e.tournament << "," << e.iteration << "," << e.generation << "," << avgFit.str() << "," << bestFit.str() << "," << e.bestGenome << "," << e.solutionFound << "," << e.solutions << "," << e.diversityMetric << endl;
    }

    file.flush();
    file.close();
}

// N,p_m,p_c,tournament_size,iteration,generation,average_fitness,best_fitness,best_genome,solution_found,num_solutions_found,diversity_metric
Experiment::Experiment(
    int N, 
    float p_m, 
    float p_c, 
    int tournament, 
    int iteration,
    int generation,
    long double averageFitness,
    pair<long double, string> bestPair,
    int solutions,
    vector <string> genomes
) {

    this->N = N;
    this->p_m = p_m;
    this->p_c = p_c;
    this->tournament = tournament;
    this->iteration = iteration;
    this->generation = generation;
    this->averageFitness = averageFitness;
    this->bestFitness = bestPair.first;
    this->bestGenome = bestPair.second;
    this->solutions = solutions;
    this->solutionFound = solutions > 0;
    this->diversityMetric = 0.0;
    this->diversityMetric = calculateDivMetric(genomes);
}

double Experiment::calculateDivMetric(vector <string>& genomes) 
{
    int N = genomes.size();
    if (N < 2) return 0.0;  // No diversity in a single-genome population

    vector<double> avgDistances(N, 0.0);
    double totalDiversity = 0.0;
    
    // Compute pairwise distances (only once per unique pair)
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) { // j > i avoids double-counting
            double dist = euclideanDistance(genomes[i], genomes[j]);
            avgDistances[i] += dist;
            avgDistances[j] += dist;
        }
    }

    // Compute the mean of pairwise distances for each genome
    for (int i = 0; i < N; i++) {
        avgDistances[i] /= (N - 1);  // Average per individual
        totalDiversity += avgDistances[i];
    }

    return totalDiversity / N;
}

void Experiment::Print() 
{
    cout << "---------------" << endl;
    cout << "N: " << N << endl;
    cout << "p_m: " << p_m << endl;
    cout << "p_c: " << p_c << endl;
    cout << "tournament: " << tournament << endl;
    cout << "iteration: " << iteration << endl;
    cout << "generation: " << generation << endl;
    cout << "averageFitness: " << averageFitness << endl;
    cout << "bestFitness: " << bestFitness << endl;
    cout << "bestGenome: " << bestGenome << endl;
    cout << "solutions: " << solutions << endl;
    cout << "solutionFound: " << solutionFound << endl;
    cout << "diversityMetric: " << diversityMetric << endl;
}