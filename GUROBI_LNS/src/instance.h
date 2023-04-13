#ifndef INSTANCE_H
#define INSTANCE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>

using namespace std;

namespace CSP {
    class Instance {

      public:

        const string file_name;
        const string input_dir;
        int cars;
        int options;
        int classes;
        int window_size;
        vector<int> max_car_in_option;
        vector<int> max_opt_block_size;
        vector<int> cars_in_class;
        vector<int> class_of_car;
        vector<vector<int> > car_opt_req;
        vector<double> pq;
        vector<double> max_pq_class;
        vector<double> min_pq_class;
        vector<double> ave_pq_class;
        vector<double> option_utilisation_class;
        vector<vector<double> > alpha_weights;
        vector<vector<double> > beta_weights;
        vector<int> opt_solution;
        std::vector<std::vector<bool>> optimal_solution_binary;
        double opt_obj;

        void read_instance();
        void display();
        void compute_alpha_beta();
        void read_optimal_solution();
        double compute_costs (const vector<int> &sequence, vector<vector<double>> &cost_car) const;
        double compute_objective_value(const vector<int> &sequence) const;

        void compute_window();
        int lcm(int a, int b);
        int gcd (int a, int b);

        bool get_optimal_value(std::uint32_t i, std::uint32_t j) const { return optimal_solution_binary[i][j]; }

        std::string get_file_name() const { return file_name; }

        // Created a new CSP instance
        explicit Instance(std::string file_name, std::string input_dir);



  };
}

#endif
