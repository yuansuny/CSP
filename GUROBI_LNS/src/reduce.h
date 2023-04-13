#ifndef REDUCE_H
#define REDUCE_H
#include "instance.h"

extern "C" {
#include "svm_predict_model.h"
#include "linear_svm_predict_model.h"
}

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <numeric>      // std::iota
#include <vector>
#include <cstring>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <iomanip>
#include <omp.h>
#include "gurobi_c++.h"

namespace CSP{
    class Reduce{
        const Instance& g;
        const float threshold_r = 0.5;
        const float threshold_c = -0.00;
        const std::string test_data_dir = "../test_data/";
        const std::string training_data_dir = "../train_data/";
        const std::string training_model_name = "train_model";
        const std::string test_data_name = "_test_data_";
        const std::string output_file_name = "_predicted_value_";
        const int sample_factor = 1000;
        const int iterations = 1000;
        int window_size = 10;
        int neighbor_size = 200;
        int sample_size = 0;
        int probability = 0;
        int num_thread = 0;
        bool restart = 1;
        double cutoff = 50;
        int m, n;
        double target_gap;
        std::vector<std::vector<int>> samples;
        std::vector<double> objs;
        std::vector<std::vector<double>> cost_car;
        double best_obj_found;
        std::vector<int> variable_order;
        std::vector<int> best_sol_found;
        std::vector<std::vector<double>> best_binary_values;
        std::vector<std::vector<bool>> predicted_value;
        std::vector<std::vector<float>> ranking_scores;
        std::vector<std::vector<float>> objective_scores;
        std::vector<std::vector<float>> corr_xy;
        std::vector<std::vector<float>> corr_xr;
        std::vector<float> min_cbm_i, min_cbm_j;
        std::vector<float> min_rcm_i, min_rcm_j;
        std::vector<float> max_rbm_i, max_rbm_j;
        std::vector<float> max_obm_i, max_obm_j;
        std::vector<float> max_cost_i, max_cost_j;
        int num_variable_left;

        GRBEnv *env = new GRBEnv();
        GRBModel model = GRBModel(*env);
        vector<vector<GRBVar> > zjt;
        vector<vector<GRBVar> > yjt;
        vector<vector<GRBVar> > xit;

        void constructing_test_data();
        void loading_output_data();
        void initialize_parameters();
        void initialize_MIP_model();
        double get_wall_clock();


//        // Prints the integer solution obtained by CPLEX to stdout.
//        void print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x) const;
        public:
            // Builds a solver for graph g.
            explicit Reduce(const Instance& g);
            void multi_thread_random_sampling();
            void multi_thread_local_sampling();
            void removing_variables_rbm();
            void removing_variables_cbm();
            void removing_variables_svm();
            void removing_variables_svm_linear();
            void repair();
            double get_objective_value_found() const { return best_obj_found; }
            bool get_predicted_value(int i, int j) const { return predicted_value[i][j]; }
            float get_cbm_value(int i, int j) const { return corr_xy[i][j]; }
            float get_rcm_value(int i, int j) const { return corr_xr[i][j]; }
            float get_rbm_value(int i, int j) const { return ranking_scores[i][j]; }
            float get_obm_value(int i, int j) const { return objective_scores[i][j]; }
            float get_min_cbm_i(int i) const { return min_cbm_i[i]; }
            float get_min_cbm_j(int j) const { return min_cbm_j[j]; }
            float get_min_rcm_i(int i) const { return min_rcm_i[i]; }
            float get_min_rcm_j(int j) const { return min_rcm_j[j]; }
            float get_max_rbm_i(int i) const { return max_rbm_i[i]; }
            float get_max_rbm_j(int j) const { return max_rbm_j[j]; }
            float get_max_obm_i(int i) const { return max_obm_i[i]; }
            float get_max_obm_j(int j) const { return max_obm_j[j]; }
            float get_max_cost_i(int i) const { return max_cost_i[i]; }
            float get_max_cost_j(int j) const { return max_cost_j[j]; }
            float get_cost_car(int i, int j) const { return cost_car[i][j]; }

            void compute_reduced_problem_size();
            void compute_correlation_based_measure();
            void compute_ranking_based_measure();
            void compute_objective_based_measure();
            void compute_ranking_correlation_measure();
            void compute_cost_measure();
            int get_reduced_problem_size() const { return num_variable_left; }
            void solve_csp_gurobi(int low, int high);
            void solve_csp_gurobi_sampling(int idx, int low, int high);
            void solver_LNS_GUROBI(double time_limit);
            void solver_LNS_sampling();

    };
}

#endif
