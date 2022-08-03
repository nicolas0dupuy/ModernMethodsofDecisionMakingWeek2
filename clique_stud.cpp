#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <cstdlib>
using namespace std;


class MaxCliqueProblem
{
public:
    static int GetRandom(int a, int b)
    {
        static mt19937 generator;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }

    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        string line;
        int vertices = 0, edges = 0;
        while (getline(fin, line))
        {
            if (line[0] == 'c')
            {
                continue;
            }

            stringstream line_input(line);
            char command;
            if (line[0] == 'p')
            {
                string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
            }
            else
            {
                int start, finish;
                line_input >> command >> start >> finish;
                // Edges in DIMACS file can be repeated, but it is not a problem for our sets
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
            }
        }
    }

    void FindClique(int randomization, int iterations)
    {
        static mt19937 generator;
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            vector<int> clique;
            vector<int> candidates(neighbour_sets.size());
            for (int i = 0; i < neighbour_sets.size(); ++i)
            {
                candidates[i] = i;
            }
            shuffle(candidates.begin(), candidates.end(), generator);
            while (! candidates.empty())
            {
                int last = candidates.size() - 1;
                int rnd = GetRandom(0, min(randomization - 1, last));
                int vertex = candidates[rnd];
                clique.push_back(vertex);
                for (int c = 0; c < candidates.size(); ++c)
                {
                    int candidate = candidates[c];
                    if (neighbour_sets[vertex].count(candidate) == 0)
                    {
                        // Move the candidate to the end and pop it
                        swap(candidates[c], candidates[last]);
                        candidates.pop_back();
                        --c;
                    }
                }
                shuffle(candidates.begin(), candidates.end(), generator);
            }
            if (clique.size() > best_clique.size())
            {
                best_clique = clique;
            }
        }
    }

    void deg_order(vector<int> & edges)
    {
        
        vector<int> ord(edges.size());
        for (int i = 0; i < edges.size(); i++)
        {
            ord[i] = neighbour_sets[edges[i]].size();
        }

        while(1)
        {
            if (ord.size() < 2)
                break;
            int is_ordered = 1;
            for (int i = 0; i < ord.size() - 1; i++)
            {
                if (ord[i] < ord[i+1])
                {
                    is_ordered = 0;
                    int tmp = ord[i];
                    ord[i] = ord[i+1];
                    ord[i+1] = tmp;
                    tmp = edges[i];
                    edges[i] = edges[i+1];
                    edges[i+1] = tmp;

                }
            }
            if (is_ordered)
                break;
        }
    }

    void MyFindClique(int iterations, int upweight_unvisited = 0, float alpha1 = 1.1, float alpha2 = 0.9, int verbose = 0, int test = 0)
        {
            /*
            upweight_unvisited : increment at the end of each iteration of the weights for unvisited vertices
            scale_iter : coef to add importance to first scores (int)(scale_iter * (iterations - current iteration)/iterations)
            */
            static mt19937 generator;
            // Variation inspired by Adaptative starting ... FIND BACK QUOTE NAME !!!
            // We define a weight vector for the probability of an edge to be picked
            vector<int> weights(neighbour_sets.size());
            for (int i=0; i < weights.size(); i++)
                // initialize the probability to the edge degree
                // weights[i] = neighbour_sets[i].size() + 1; 
                weights[i] = 1;

            for (int iteration = 0; iteration < iterations; ++iteration)
            {
                vector<int> clique;
                vector<int> candidates(neighbour_sets.size());
                int last = candidates.size() - 1;
                for (int i = 0; i < neighbour_sets.size(); ++i)
                {
                    candidates[i] = i;
                }
                // Removing edges of degree < best_clique size (sure not to belong to the best clique)
                for (int i = 0 ; i < candidates.size(); i++)
                {
                    if (neighbour_sets[i].size() < best_clique.size())
                    {
                        int last = candidates.size() - 1;
                        swap(candidates[i], candidates[last]);
                        candidates.pop_back();
                        --i;
                    }
                }
                
                while (! candidates.empty())
                {
                    // Get the total weight for remaining candidates
                    int total_weights = 0;
                    last = candidates.size() - 1;
                    for (int i = 0; i < candidates.size(); i++)
                        total_weights += weights[candidates[i]];
                    int t = GetRandom(0, total_weights);
                    int sum = 0;
                    int rnd = 0;
                    // Find the last edge so that the sum of weights before it is < to t
                    // Leads to a weighted probability
                    shuffle(candidates.begin(), candidates.end(), generator);
                    sum += weights[candidates[rnd]];
                    while (sum < t)
                    {
                        rnd += 1;
                        sum += weights[candidates[rnd]];
                    }
                    int vertex = candidates[rnd];
                    clique.push_back(vertex);
                    for (int c = 0; c < candidates.size(); ++c)
                    {
                        int candidate = candidates[c];
                        if (neighbour_sets[vertex].count(candidate) == 0)
                        {
                            // Move the candidate to the end and pop it
                            swap(candidates[c], candidates[last]);
                            candidates.pop_back();
                            --c;
                        }
                    }
                }
                if (verbose)
                {
                    cout << "Current clique : " ;
                    for (int i : clique)
                        cout << i << " ";
                    cout << endl;
                }
                /*
                if (test)
                {
                    best_clique = clique;

                    if (! Check())
                    {
                        cout << "Not a clique ! Stopping" << endl;
                        return;
                    }
                }
                */
                if (clique.size() > best_clique.size())
                {
                    best_clique = clique;
                }
                // increment weights of unvisited edges (all weights updated, then the next loop downgrade visited edges)
                // allows for exploration
                int score = (clique.size() - best_clique.size());
                if (score > 0)
                {
                    score *= alpha1;
                    alpha1 *= alpha1;
                }
                if (score < 0)
                {
                    score *= alpha2;
                    alpha2 *= alpha2;
                }
                for (int i = 0; i < weights.size(); i++)
                {
                    int is_in_clique = 0;
                    for (int j = 0; j < clique.size(); j++)
                    {
                        if (i == clique[j])
                            {
                                is_in_clique = 1;
                                weights[clique[i]] += score;
                            }
                    }
                    if (! is_in_clique)
                        weights[i] += upweight_unvisited;

                }
            }
        }


    const vector<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        if (unique(best_clique.begin(), best_clique.end()) != best_clique.end())
        {
            cout << "Duplicated vertices in the clique\n";
            return false;
        }
        for (int i : best_clique)
        {
            for (int j : best_clique)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not a clique\n";
                    return false;
                }
            }
        }
        return true;
    }

private:
    vector<unordered_set<int>> neighbour_sets;
    vector<int> best_clique;
};

int main(int argc, char *argv[])
{
    int iterations = 100;
    int upweight_unvisited = 1;
    int scale_iter = 1;
    int max_weight = 1000;
    int verbose = 0;
    int test = 0;
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            if(i == 1)
                iterations = atoi(argv[i]);
            if(i == 2)
                upweight_unvisited = atoi(argv[i]);
            if(i == 3)
                scale_iter = atoi(argv[i]);
            if (i == 4)
                max_weight = atoi(argv[i]);
            if (i == 5)
                verbose = atoi(argv[i]);
            if (i == 6)
                test = atoi(argv[i]);
            }
    }
    cout << "USAGE :" << endl;
    cout << "1st argument : Niter : " << iterations << endl;
    cout << "2nd argument : increment value for unvisited edges : " << upweight_unvisited<< endl;
    cout << "3rd argument : scale_iter : " << scale_iter << endl;
    cout << "3rd argument : max_weight : " << max_weight << endl;
    cout << "4th argument : verbose" << endl;
    cout << "5th argument : test (unimplemented yet)" << endl;
    cout << "    ****    ****    ****    ****    ****    ****    ****" << endl;
    /*
    cout << "Number of iterations: ";
    cin >> iterations;
    */
    /*
    int randomization;
    cout << "Randomization: ";
    cin >> randomization;
    */
    /*
    vector<string> files = { "C125.9.clq", "johnson8-2-4.clq", "johnson16-2-4.clq", "MANN_a9.clq", "MANN_a27.clq",
        "p_hat1000-1.clq", "keller4.clq", "hamming8-4.clq", "brock200_1.clq", "brock200_2.clq", "brock200_3.clq", "brock200_4.clq",
        "gen200_p0.9_44.clq", "gen200_p0.9_55.clq", "brock400_1.clq", "brock400_2.clq", "brock400_3.clq", "brock400_4.clq",
        "MANN_a45.clq", "sanr400_0.7.clq", "p_hat1000-2.clq", "p_hat500-3.clq", "p_hat1500-1.clq", "p_hat300-3.clq", "san1000.clq",
        "sanr200_0.9.clq" };
    */
    vector<string> files = {"brock200_1.clq"};
    ofstream fout("clique.csv");
    fout << "File; Clique; Time (sec)\n";
    for (string file : files)
    {
        MaxCliqueProblem problem;
        problem.ReadGraphFile(file);
        clock_t start = clock();
        problem.MyFindClique(iterations, upweight_unvisited, scale_iter, max_weight, verbose, test);
        if (! problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }
        fout << file << "; " << problem.GetClique().size() << "; " << double(clock() - start) / 1000 << '\n';
        cout << file << ", result - " << problem.GetClique().size() << ", time - " << double(clock() - start) / 1000 << '\n';
        cout << "{" ;
        for (int i : problem.GetClique())
            cout << i << " ";
        cout << "}" << endl;
    }
    fout.close();
    return 0;
}
