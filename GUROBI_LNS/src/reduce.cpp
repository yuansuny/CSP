#include "reduce.h"
#include <cmath>
#include <limits>

namespace CSP {
    Reduce::Reduce(const Instance& g) : g{g}{
        initialize_parameters();
    }

    void Reduce::solver_LNS_sampling(){
        std::string output_sol_filename = "../results/" + g.file_name + "_res_sol.csv";
        std::ofstream output_file_sol (output_sol_filename, std::ios::trunc);
        if (output_file_sol.is_open()){
            output_file_sol << "Sample solution" << "\n";
        } else{
            std::cout << "Cannot open the output file " + output_sol_filename << "\n";
            return;
        }
        output_file_sol.close();

        double w0 = get_wall_clock();
        initialize_MIP_model();
        multi_thread_random_sampling();

        int high;
        double best_objective_pre;

        for (int i = 0; i < sample_size; ++i){
            std::cout << "sample " << i << std::endl;
            best_objective_pre = 100000000;
            while((best_objective_pre - objs[i])/objs[i] > 0.00005){
                best_objective_pre = objs[i];
                for (int j = 0; j < n; j = j + window_size/2){
                    if (j + window_size < n){
                        high = j + window_size;
                    }else{
                        high = n;
                    }
                    solve_csp_gurobi_sampling(i, j, high);
                    if (high == n){
                        break;
                    }
                }
            }

            if (objs[i] < best_obj_found){
                for (int j = 0; j < n; ++j){
                    best_sol_found[j] = samples[i][j];
                }
                best_obj_found = objs[i];
            }
            std::cout << "best objective value found so far is " << best_obj_found << std::endl;
            output_file_sol.open(output_sol_filename, std::ios::app);
            for (int j = 0; j < n; ++j){
                output_file_sol << samples[i][j] << ", ";
            }
            output_file_sol << "\n";
            output_file_sol.close();
        }
    }

    void Reduce::solver_LNS_GUROBI(double time_limit) {

        double w0 = get_wall_clock();
        initialize_MIP_model();
        multi_thread_random_sampling();

        model.set(GRB_DoubleParam_TimeLimit, cutoff);
        solve_csp_gurobi(0, n);


        std::random_device rd;
        std::mt19937 gen(rd());

        int high;
        int iter = 0;
        double best_obj_pre;
        while(get_wall_clock() - w0 < time_limit){
            best_obj_pre = best_obj_found;
            std::shuffle(variable_order.begin(), variable_order.end(), gen);
            for (int i = 0; i < n; i = i + window_size/2){
                if (i + window_size < n){
                    high = i + window_size;
                }else{
                    high = n;
                }
                if (get_wall_clock() - w0 < time_limit){
                    model.set(GRB_DoubleParam_TimeLimit, time_limit - get_wall_clock() + w0);
                    solve_csp_gurobi(i, high);
                }
                if (high == n){
                    break;
                }
            }
            std::cout << "Iteration: " << ++iter << ", time used is: " << get_wall_clock() - w0 << \
                          ", best objective value found is: " << best_obj_found << std::endl;

            if ((best_obj_pre - best_obj_found)/best_obj_found < 0.005){
                window_size = window_size + 1;
                std::cout << "window size is : " <<  window_size << std::endl;
            }

            if (window_size > n){
                return;
            }
        }
    }

//        for (int i = 0; i < iterations; ++i){
//            std::cout << "iteration : " << i+1 << std::endl;
//            if (restart == 1){
//                multi_thread_random_sampling();
//            }else{
//                multi_thread_local_sampling();
//            }
//
//            removing_variables_svm();
//
////            model.getEnv().set(GRB_DoubleParam_MIPGap, 0.005);
////            target_gap = 0.02*iterations - 0.02*i;
//            solve_csp_gurobi();
//
////            if (i == iterations - 1){
////                for (int w = 0; w < m; ++w){
////                    for (int j = 0; j < n; ++j){
////                        predicted_value[w][j] = 1;
////                    }
////                }
////            }
//
//            // solve the reduced problem by Gurobi
////            GRBLinExpr rtot = 0;
////            for (int w = 1; w < m + 1; ++w){
////                for (int t = 1; t < n + 1; ++t){
////                    if (predicted_value[w-1][t-1] < 0.5){
////                        rtot += xit[w][t];
////                    }
////                }
////            }
////            for (int j = 0; j < 2; ++j){
////                if (j == 0){
////                    model.addConstr(rtot <= 5, "fixing_constraint");
////                } else {
////                    model.remove(model.getConstrByName("fixing_constraint"));
////                    if (j < 2){
////                        model.addConstr(rtot >= 6, "fixing_constraint");
////                    }
//////                    else{
//////                        model.addConstr(rtot >= j, "fixing_constraint");
//////                    }
////                }
////                model.update();

    void Reduce::solve_csp_gurobi_sampling(int idx, int low, int high){
        for (int t = 1; t < n + 1; ++t){
            if (t-1<low || t-1>=high){
                for (int w = 1; w < m + 1; ++w){
                    if (w == samples[idx][t-1] + 1){
                        xit[w][t].set(GRB_DoubleAttr_LB, 1);
                    }else{
                        xit[w][t].set(GRB_DoubleAttr_UB, 0);
                    }
                }
            }
        }
        model.update();

//      warm start MIP
        for (int t = 1; t < n + 1; ++t){
            for (int w = 1; w < m + 1; ++w){
                if (w == samples[idx][t-1] + 1){
                    xit[w][t].set(GRB_DoubleAttr_Start, 1);
                }else{
                    xit[w][t].set(GRB_DoubleAttr_Start, 0);
                }
            }
        }
        model.update();

        model.optimize();
        double objVal = model.get(GRB_DoubleAttr_ObjVal);
        double objBound = model.get(GRB_DoubleAttr_ObjBound);

        // copy the best found solution
        for (int t = 1; t < n+1; ++t) {
            for (int w = 1; w < m+1; ++w){
                if (xit[w][t].get(GRB_DoubleAttr_X) > 0.5){
                    samples[idx][t-1] = w - 1;
                }
            }
        }
        objs[idx] = objVal;

        // adding removed varaibles back to the model
        for (int t = 1; t < n + 1; ++t){
            for (int w = 1; w < m + 1; ++w){
                xit[w][t].set(GRB_DoubleAttr_UB, 1);
                xit[w][t].set(GRB_DoubleAttr_LB, 0);
            }
        }
        model.update();
        model.reset();
    }


    void Reduce::solve_csp_gurobi(int low, int high){

//        std::cout << "low is : " << low << " , high is : " << high << std::endl;

        // removing variables by setting upper bound = lower bound = 0
        for (int t = 1; t < n + 1; ++t){
            for (int w = 1; w < m + 1; ++w){
                if (t-1<low || t-1>=high){
//                    if (best_binary_values[w-1][t-1] < 0.5){
//                        xit[w][t].set(GRB_DoubleAttr_UB, 0);
//                    }else{
//                        xit[w][t].set(GRB_DoubleAttr_LB, 1);
//                    }
                    if (best_binary_values[w-1][variable_order[t-1]] < 0.5){
                        xit[w][variable_order[t-1]+1].set(GRB_DoubleAttr_UB, 0);
                    }else{
                        xit[w][variable_order[t-1]+1].set(GRB_DoubleAttr_LB, 1);
                    }
                }
            }
        }
        model.update();

//         warm start MIP
        for (int w = 1; w < m + 1; ++w){
            for (int t = 1; t < n + 1; ++t){
                xit[w][t].set(GRB_DoubleAttr_Start, best_binary_values[w-1][t-1]);
            }
        }

        model.optimize();
        double objVal = model.get(GRB_DoubleAttr_ObjVal);
        double objBound = model.get(GRB_DoubleAttr_ObjBound);

        // copy the best found solution
        if (objVal < best_obj_found){
            for (int t = 1; t < n+1; ++t) {
                for (int w = 1; w < m+1; ++w){
                    best_binary_values[w-1][t-1] = xit[w][t].get(GRB_DoubleAttr_X);
                    if (xit[w][t].get(GRB_DoubleAttr_X) > 0.5){
                        best_sol_found[t-1] = w - 1;
                    }
                }
            }
            best_obj_found = objVal;
        }

        // adding removed varaibles back to the model
        for (int t = 1; t < n + 1; ++t){
            for (int w = 1; w < m + 1; ++w){
                xit[w][t].set(GRB_DoubleAttr_UB, 1);
                xit[w][t].set(GRB_DoubleAttr_LB, 0);
            }
        }
        model.update();
        model.reset();

//        std::cout << "best objective found by Gurobi is : " << objVal << std::endl;

    }

    void Reduce::initialize_MIP_model(){
//        GRBEnv *env;
//        env = new GRBEnv();
//        model = GRBModel(*env);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_StringAttr_ModelName, "CSP");
        // Create variables
        zjt.resize(g.options+1);
        yjt.resize(g.options+1);
        xit.resize(m+1);
        for (int w = 1; w < g.options+1; w++) {
            zjt[w].resize(n+1);
            yjt[w].resize(n+1);
            for(int t = 0; t < n+1; t++){
                zjt[w][t] = model.addVar(0,n,0,GRB_INTEGER);
                yjt[w][t] = model.addVar(0,n,0,GRB_INTEGER);
            }
        }
        for (int w = 1; w < m + 1; w++) {
            xit[w].resize(n+1);
            for(int t = 1; t < n+1; t++){
                xit[w][t] = model.addVar(0,1,0,GRB_BINARY);
            }
        }
        model.update();

        // only one car sequenced at one point
        for (int i = 1; i < n+1; i++) {
            GRBLinExpr rtot = 0;
            for(int t = 1; t < m+1; t++){
                rtot += xit[t][i];
            }
            model.addConstr(rtot == 1, "");
        }
        model.update();

        // the number of an option used across the time points = number of cars that need that option
        for (int i = 1; i < m+1; i++) {
            GRBLinExpr rtot = 0;
            for(int t = 1; t < n+1; t++){
                rtot += xit[i][t];
            }
            model.addConstr(rtot == g.cars_in_class[i-1], "");
        }
        model.update();

        // map between zjt, yjt and xit
        for (int t = 1; t < n+1; t++) {
            for (int j = 1; j < g.options+1; j++) {
                GRBLinExpr rtot = 0;
                int min_p = 1;
                if(g.max_opt_block_size[j-1] < t) min_p = t - g.max_opt_block_size[j-1] + 1;
                for(int i = 1; i < m+1; i++){
                    for(int k = min_p; k < t+1; k++) {
                        rtot += g.car_opt_req[i-1][j-1] * xit[i][k];
                    }
                }
                if(t < g.max_car_in_option[j-1])
                    model.addConstr(zjt[j][t] - yjt[j][t] + rtot == t, "");
                else
                    model.addConstr(zjt[j][t] - yjt[j][t] + rtot == g.max_car_in_option[j-1], "");
            }
        }
        model.update();

        // set time limit
//        model.set(GRB_DoubleParam_TimeLimit, cutoff);
        // set number of thread
        model.set(GRB_IntParam_Threads, omp_get_max_threads());
        // set random seed
        srand (time(NULL));
        model.set(GRB_IntParam_Seed, rand()%10000);
        // MIP focus
//        model.set(GRB_IntParam_MIPFocus, 1);

//      model.set(GRB_IntParam_BarIterLimit, 30);
      model.set(GRB_DoubleParam_MIPGap, 0.005);
      model.set(GRB_IntParam_Presolve, 1);
//      model.set(GRB_IntParam_Threads, 4);
//      model.set(GRB_IntParam_BarOrder, 1);
//      model.set(GRB_IntParam_Cuts,-1); // aggressive cut generation
//      model.set(GRB_IntParam_CutPasses,200); // maximum number of cutting plane passes for the root cut generation
//      // Not sure what the below mean.
//      model.set(GRB_DoubleParam_BarConvTol,0.0001);
//      model.set(GRB_IntParam_PreSparsify,1);
//      model.set(GRB_IntParam_PrePasses,2);


        // the objective
        // find the objective
        GRBLinExpr tot=0;
        for(int i = 1; i < g.options+1; i++){
            for(int t = 1; t < n+1; t++){
                tot += g.alpha_weights[t-1][i-1] * yjt[i][t] + g.beta_weights[t-1][i-1] * zjt[i][t];
            }
        }
        //set the objective
        model.setObjective(tot,GRB_MINIMIZE);
        model.update();
    }

    void Reduce::initialize_parameters() {
        n = g.cars;
        m = g.classes;
//        window_size = g.window_size;
        window_size = 30;
        neighbor_size = n;
        variable_order = std::vector<int>(n);
        std::iota(variable_order.begin(), variable_order.end(), 0);
//        sample_size = sample_factor * m;
        sample_size = 100;
        std::cout << "sample size is " << sample_size << "\n";
        num_thread = omp_get_max_threads();
        std::cout << "threads used: " << num_thread << "\n";
        cutoff = (double) n / num_thread;
//        cutoff = 10;
//        cutoff = (double) n;
        samples = std::vector<std::vector<int>>(sample_size, std::vector<int>(n));
        best_sol_found = std::vector<int>(n);
        best_obj_found = 9999999;
        objs = std::vector<double>(sample_size, 0.0);
        predicted_value = std::vector<std::vector<bool>>(m, std::vector<bool>(n, 0));
        best_binary_values = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));

        // cost score
        cost_car = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));
        max_cost_i = std::vector<float>(m, 0.0);
        max_cost_j = std::vector<float>(n, 0.0);

        // ranking-based measure
        ranking_scores = std::vector<std::vector<float>>(m, std::vector<float>(n, 0.0));
        max_rbm_i = std::vector<float>(m, 0.0);
        max_rbm_j = std::vector<float>(n, 0.0);

        // objective-based measure
        objective_scores = std::vector<std::vector<float>>(m, std::vector<float>(n, 0.0));
        max_obm_i = std::vector<float>(m, 0.0);
        max_obm_j = std::vector<float>(n, 0.0);

        // correlation-based measure
        corr_xy = std::vector<std::vector<float>>(m, std::vector<float>(n, 0.0));
        min_cbm_i = std::vector<float>(m, 1.0);
        min_cbm_j = std::vector<float>(n, 1.0);

        // ranking correlation-based measure
        corr_xr = std::vector<std::vector<float>>(m, std::vector<float>(n, 0.0));
        min_rcm_i = std::vector<float>(m, 1.0);
        min_rcm_j = std::vector<float>(n, 1.0);
    }

//    void Reduce::multi_thread_local_sampling() {
//        for (int i = 0; i < m; ++i){
//            for (int j = 0; j < n; ++j){
//                cost_car[i][j] = 0.0;
//            }
//        }
//        #pragma omp parallel
//        {
//            int threadnum = omp_get_thread_num();
//            int low = sample_size * threadnum / num_thread;
//            int high = sample_size * (threadnum+1) / num_thread;
//            std::vector<int> best_sol_local = std::vector<int>(n);
//            std::random_device rd;
//            std::mt19937 gen(rd()+threadnum);
//            std::uniform_int_distribution<> dist(0, n-1);
//            srand (time(NULL));
//            int idx1, idx2;
//            double val, obj_val_local;
////            for (int i = 0; i < n; ++i){
////                best_sol_local[i] = best_sol_sampling[i];
////            }
//
//            for(int i = low; i < high; ++i) {
//                if (i == low && threadnum == 0){
//                    for (int j = 0; j < n; ++j){
//                        samples[i][j] = best_sol_sampling[j];
//                    }
////                    objs[i] = best_obj_sampling;
//                }else{
////                for (int j = 0; j < n; ++j){
////                    best_sol_local[j] = samples[i][j];
////                }
////                if (i > low){
//                    std::shuffle(samples[i].begin(), samples[i].end(), gen);
////                    for (int j = 0; j < neighbor_size; ++j){
////                        idx1 = rand() % n;
////                        idx2 = rand() % n;
////                        val = best_sol_local[idx1];
////                        best_sol_local[idx1] = best_sol_local[idx2];
////                        best_sol_local[idx2] = val;
////                    }
//                }
//                objs[i] = g.compute_objective_value(samples[i]);
////                obj_val_local = g.compute_objective_value(best_sol_local);
////                if (obj_val_local < objs[i]){
////                    for (int j = 0; j < n; ++j){
////                        samples[i][j] = best_sol_local[j];
////                    }
////                }
//            }
//        }
//
//        int idx = 0;
//        double min_obj = objs[0];
//        for (int i = 1; i < sample_size; ++i){
//            if (objs[i] < min_obj){
//                min_obj = objs[i];
//                idx = i;
//            }
//        }
//        std::cout << "best obj in local sampling is " << objs[idx] << "\n";
//        for (int j = 0; j < n; ++j){
//            best_sol_sampling[j] = samples[idx][j];
//        }
//        best_obj_sampling = objs[idx];
//    }

    void Reduce::multi_thread_local_sampling() {
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
            int low = sample_size * threadnum / num_thread;
            int high = sample_size * (threadnum+1) / num_thread;
            std::vector<int> best_sol_local = std::vector<int>(n);
            double best_obj_local;
            std::random_device rd;
            std::mt19937 gen(rd()+threadnum);
            std::uniform_int_distribution<> dist1(low, high-1);
            std::uniform_int_distribution<> dist2(0, n-1);
//            srand(threadnum);
            int idx1, idx2, idx;
            double val;
            best_obj_local = best_obj_found;
            for (int j = 0; j < n; ++j){
                best_sol_local[j] = best_sol_found[j];
            }
//            for(int i = 0; i < m / num_thread; ++i) {
//            for(int i = 0; i < 1; ++i) {
            for(int i = low; i < high; ++i) {
//                idx = dist1(gen);
                idx = i;
                for (int j = 0; j < n; ++j){
                    samples[idx][j] = best_sol_local[j];
                }
                for (int j = 0; j < neighbor_size; ++j){
                    idx1 = dist2(gen);
                    idx2 = dist2(gen);
//                    idx1 = rand()%n;
//                    idx2 = rand()%n;
//
//                    val = samples[idx][idx1];
//                    samples[idx][idx1] = samples[idx][idx2];
//                    samples[idx][idx2] = val;
                    std::swap(samples[idx][idx1], samples[idx][idx2]);
                }
                objs[idx] = g.compute_objective_value(samples[idx]);
                if (objs[idx] < best_obj_local){
                    best_obj_local = objs[idx];
                    for (int j = 0; j < n; ++j){
                        best_sol_local[j] = samples[idx][j];
                    }
                }
            }
            #pragma omp critical
            if(best_obj_local < best_obj_found){
                best_obj_found = best_obj_local;
                for (int j = 0; j < n; ++j){
                    best_sol_found[j] = best_sol_local[j];
                }
            }
        }
    }

//    void Reduce::multi_thread_local_sampling() {
//        #pragma omp parallel
//        {
//            int threadnum = omp_get_thread_num();
//            int low = sample_size * threadnum / num_thread;
//            int high = sample_size * (threadnum+1) / num_thread;
//            std::vector<int> sol_local = std::vector<int>(n);
//            std::vector<int> order = std::vector<int>(n);
//            std::iota(order.begin(), order.end(), 0);
//            std::random_device rd;
//            std::mt19937 gen(rd()+threadnum);
//            srand (time(NULL));
//            int idx1, idx2;
//            double obj_val_local;
//
//            for(int i = low; i < high; ++i) {
//                for (int j = 0; j < n; ++j){
//                    sol_local[j] = samples[i][j];
//                }
//                std::shuffle(order.begin(), order.end(), gen);
//                for (int j = 0; j < n; ++j){
//                    idx1 = j;
//                    idx2 = order[j];
//                    if (sol_local[idx1] != sol_local[idx2]){
//                        std::swap(sol_local[idx1], sol_local[idx2]);
//                        obj_val_local = g.compute_objective_value(sol_local);
//                        if (obj_val_local < objs[i]){
//                            objs[i] = obj_val_local;
//                        }else{
//                            std::swap(sol_local[idx1], sol_local[idx2]);
//                        }
//                    }
//                }
//                for (int j = 0; j < n; ++j){
//                    samples[i][j] = sol_local[j];
//                }
//            }
//        }
//
//        int idx = 0;
//        double min_obj = objs[0];
//        for (int i = 1; i < sample_size; ++i){
//            if (objs[i] < min_obj){
//                min_obj = objs[i];
//                idx = i;
//            }
//        }
//        std::cout << "best obj in local sampling is " << objs[idx] << "\n";
//        for (int j = 0; j < n; ++j){
//            best_sol_sampling[j] = samples[idx][j];
//        }
//        best_obj_sampling = objs[idx];
//    }


//    void Reduce::multi_thread_local_sampling() {
//        #pragma omp parallel
//        {
//            int threadnum = omp_get_thread_num();
//            int low = sample_size * threadnum / num_thread;
//            int high = sample_size * (threadnum+1) / num_thread;
//            std::vector<int> sol_local = std::vector<int>(n);
//            std::vector<int> order = std::vector<int>(n);
//            std::iota(order.begin(), order.end(), 0);
//            std::random_device rd;
//            std::mt19937 gen(rd()+threadnum);
//            srand (time(NULL));
//            int idx1, idx2;
//            double obj_val_local;
//
//            for(int i = low; i < high; ++i) {
//                for (int j = 0; j < n; ++j){
//                    sol_local[j] = samples[i][j];
//                }
//                std::shuffle(order.begin(), order.end(), gen);
//                for (int j = 0; j < n; ++j){
//                    for (int k = j+1; k < n; ++k){
//                        idx1 = order[j];
//                        idx2 = order[k];
//                        if (sol_local[idx1] != sol_local[idx2]){
//                            std::swap(sol_local[idx1], sol_local[idx2]);
//                            obj_val_local = g.compute_objective_value(sol_local);
//                            if (obj_val_local < objs[i]){
//                                objs[i] = obj_val_local;
//                            }else{
//                                std::swap(sol_local[idx1], sol_local[idx2]);
//                            }
//                        }
//                    }
//                }
//                for (int j = 0; j < n; ++j){
//                    samples[i][j] = sol_local[j];
//                }
//            }
//        }
//
//        int idx = 0;
//        double min_obj = objs[0];
//        for (int i = 1; i < sample_size; ++i){
//            if (objs[i] < min_obj){
//                min_obj = objs[i];
//                idx = i;
//            }
//        }
//        std::cout << "best obj in local sampling is " << objs[idx] << "\n";
//        for (int j = 0; j < n; ++j){
//            best_sol_sampling[j] = samples[idx][j];
//        }
//        best_obj_sampling = objs[idx];
//    }

    void Reduce::multi_thread_random_sampling() {
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
            int low = sample_size * threadnum / num_thread;
            int high = sample_size * (threadnum+1) / num_thread;
            std::random_device rd;
            std::mt19937 gen(rd()+threadnum);
            srand (time(NULL));
            for(int i = low; i < high; ++i) {
                for (int j = 0; j < n; ++j){
                    samples[i][j] = g.class_of_car[j];
                }
                std::shuffle(samples[i].begin(), samples[i].end(), gen);
                objs[i] = g.compute_objective_value(samples[i]);
            }
        }

        int idx = -1;
        for (int i = 0; i < sample_size; ++i){
            if (objs[i] < best_obj_found){
                best_obj_found = objs[i];
                idx = i;
            }
        }
        if (idx != -1){
            std::cout << "best obj in random sampling is " << objs[idx] << "\n";
            for (int j = 0; j < n; ++j){
                best_sol_found[j] = samples[idx][j];
                best_binary_values[best_sol_found[j]][j] = 1;
            }
        }
    }


    void Reduce::compute_cost_measure() {
        std::cout << "computing objective measure " << std::endl;
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                cost_car[i][j] = 0.0;
            }
        }
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
            int low = sample_size * threadnum / num_thread;
            int high = sample_size * (threadnum+1) / num_thread;
            std::vector<std::vector<double>> cost_car_local = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));
            for(int i = low; i < high; ++i) {
                objs[i] = g.compute_costs(samples[i], cost_car_local);
            }
            #pragma omp critical
            for (int i = 0; i < m; ++i){
                for (int j = 0; j < n; ++j){
                    cost_car[i][j] += cost_car_local[i][j] / sample_size;
                }
            }
        }

        for (int i = 0; i < m; ++i){
            max_cost_i[i] = 0.0;
        }
        for (int j = 0; j < n; ++j){
            max_cost_j[j] = 0.0;
        }
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                if (cost_car[i][j] > max_cost_i[i]){
                    max_cost_i[i] = cost_car[i][j];
                }
                if (cost_car[i][j] > max_cost_j[j]){
                    max_cost_j[j] = cost_car[i][j];
                }
            }
        }
    }


    void Reduce::removing_variables_rbm() {
        std::cout << "problem reduction using ranking-based measure with threshold " << threshold_r << std::endl;
        compute_ranking_based_measure();
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                if (ranking_scores[i][j] < threshold_r){
                    predicted_value[i][j] = 0;
                } else{
                    predicted_value[i][j] = 1;
                }
            }
        }
        repair();
        compute_reduced_problem_size();
    }

    void Reduce::compute_ranking_based_measure() {
        std::cout << "computing ranking based measure " << std::endl;
        std::vector<int> sort_idx(sample_size);
        std::iota(sort_idx.begin(), sort_idx.end(), 0);
        std::vector<double> objs_copy(objs);
        std::sort(sort_idx.begin(), sort_idx.end(), [&objs_copy](int i1, int i2) {return objs_copy[i1] < objs_copy[i2];});

        std::vector<int> rank = std::vector<int>(sample_size);
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                ranking_scores[i][j] = 0.0;
            }
        }
        for (int i = 0; i < sample_size; i++){
            rank[sort_idx[i]] = i;
        }

        #pragma omp parallel
        {
//            std::cout << "computing ranking scores from thread " << omp_get_thread_num() << std::endl;
            std::vector<std::vector<float>> ranking_scores_local = std::vector<std::vector<float>>(m, std::vector<float>(n, 0.0));
            int threadnum = omp_get_thread_num();
            int low = sample_size * threadnum / num_thread;
            int high = sample_size * (threadnum+1) / num_thread;
            for (int i = low; i < high; ++i){
                for (int j = 0; j < n; ++j){
                    ranking_scores_local[samples[i][j]][j] += 1.0/(rank[i]+1);
                }
            }
            #pragma omp critical
            for (int i = 0; i < m; ++i){
                for (int j = 0; j < n; ++j){
                    ranking_scores[i][j] += ranking_scores_local[i][j];
                }
            }
        }

        for (int i = 0; i < m; ++i){
            max_rbm_i[i] = 0.0;
        }
        for (int j = 0; j < n; ++j){
            max_rbm_j[j] = 0.0;
        }
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                if (ranking_scores[i][j] > max_rbm_i[i]){
                    max_rbm_i[i] = ranking_scores[i][j];
                }
                if (ranking_scores[i][j] > max_rbm_j[j]){
                    max_rbm_j[j] = ranking_scores[i][j];
                }
            }
        }
    }


    void Reduce::compute_objective_based_measure() {
        std::cout << "computing objective based measure " << std::endl;
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                objective_scores[i][j] = 0.0;
            }
        }
        #pragma omp parallel
        {
            std::vector<std::vector<float>> objective_scores_local = std::vector<std::vector<float>>(m, std::vector<float>(n, 0.0));
            int threadnum = omp_get_thread_num();
            int low = sample_size * threadnum / num_thread;
            int high = sample_size * (threadnum+1) / num_thread;
            for (int i = low; i < high; ++i){
                for (int j = 0; j < n; ++j){
                    objective_scores_local[samples[i][j]][j] += 1.0/(objs[i]);
                }
            }
            #pragma omp critical
            for (int i = 0; i < m; ++i){
                for (int j = 0; j < n; ++j){
                    objective_scores[i][j] += objective_scores_local[i][j];
                }
            }
        }

        for (int i = 0; i < m; ++i){
            max_obm_i[i] = 0.0;
        }
        for (int j = 0; j < n; ++j){
            max_obm_j[j] = 0.0;
        }
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                if (objective_scores[i][j] > max_obm_i[i]){
                    max_obm_i[i] = objective_scores[i][j];
                }
                if (objective_scores[i][j] > max_obm_j[j]){
                    max_obm_j[j] = objective_scores[i][j];
                }
            }
        }
    }


    void Reduce::removing_variables_cbm(){
        std::cout << "problem reduction using correlation-based measure with threshold " << threshold_c << std::endl;
        compute_correlation_based_measure();
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                if (corr_xy[i][j] < threshold_c){
                    predicted_value[i][j] = 1;
                } else{
                    predicted_value[i][j] = 0;
                }
            }
        }
        repair();
        compute_reduced_problem_size();
    }

    void Reduce::compute_correlation_based_measure(){
        std::cout << "computing correlation based measure " << std::endl;
        double mean_y = 0.0;
        for (int i = 0; i < sample_size; ++i){
            mean_y += objs[i];
        }
        mean_y = mean_y/sample_size;
        std::vector<double> diff_y = std::vector<double>(sample_size);
        double variance_y = 0.0, sum_diff_y = 0.0;
        for (int i = 0; i < sample_size; ++i){
            diff_y[i] = objs[i] - mean_y;
            variance_y += diff_y[i]*diff_y[i];
            sum_diff_y += diff_y[i];
        }

        std::vector<std::vector<double>> mean_x = std::vector<std::vector<double>> (m, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> S1 = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));

        #pragma omp parallel
        {
//            std::cout << "computing correlation measure from thread " << omp_get_thread_num() << std::endl;
            std::vector<std::vector<double>> mean_x_local = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));
            std::vector<std::vector<double>> S1_local = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));

            int threadnum = omp_get_thread_num();
            int low = sample_size * threadnum / num_thread;
            int high = sample_size * (threadnum+1) / num_thread;
            double ratio = 1.0/sample_size;

            for (int i = low; i < high; ++i){
                for (int j = 0; j < n; ++j){
                    mean_x_local[samples[i][j]][j] += ratio;
                    S1_local[samples[i][j]][j] += diff_y[i];
                }
            }

            #pragma omp critical
            for (int i = 0; i < m; ++i){
                for (int j = 0; j < n; ++j){
                    mean_x[i][j] += mean_x_local[i][j];
                    S1[i][j] += S1_local[i][j];
                }
            }
        }

        std::vector<std::vector<double>> variance_x(m, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> variance_xy(m, std::vector<double>(n, 0.0));
        for (int i = 0; i < m; ++i){
            min_cbm_i[i] = 1.0;
        }
        for (int j = 0; j < n; ++j){
            min_cbm_j[j] = 1.0;
        }
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                variance_x[i][j] = mean_x[i][j]*(1-mean_x[i][j])*sample_size;
                variance_xy[i][j] = (1-mean_x[i][j])*S1[i][j] - mean_x[i][j]*(sum_diff_y - S1[i][j]);
                if (variance_x[i][j] > 0){
                    corr_xy[i][j] = variance_xy[i][j]/sqrt(variance_x[i][j]*variance_y);
                }else if (mean_x[i][j] == 0) {
                    corr_xy[i][j] = 1.0;
                } else{
                    corr_xy[i][j] = -1.0;
                }
                if (corr_xy[i][j] < min_cbm_i[i]){
                    min_cbm_i[i] = corr_xy[i][j];
                }
                if (corr_xy[i][j] < min_cbm_j[j]){
                    min_cbm_j[j] = corr_xy[i][j];
                }
            }
        }
    }

    void Reduce::compute_ranking_correlation_measure(){
        std::cout << "computing rank correlation measure " << std::endl;

        std::vector<int> sort_idx(sample_size);
        std::iota(sort_idx.begin(), sort_idx.end(), 0);
        std::vector<double> objs_copy(objs);
        std::sort(sort_idx.begin(), sort_idx.end(), [&objs_copy](int i1, int i2) {return objs_copy[i1] < objs_copy[i2];});
        std::vector<int> rank = std::vector<int>(sample_size);
        for (int i = 0; i < sample_size; i++){
            rank[sort_idx[i]] = i;
        }

        double mean_y = 0.0;
        for (int i = 0; i < sample_size; ++i){
            mean_y += rank[i];
        }
        mean_y = mean_y/sample_size;
        std::vector<double> diff_y = std::vector<double>(sample_size);
        double variance_y = 0.0, sum_diff_y = 0.0;
        for (int i = 0; i < sample_size; ++i){
            diff_y[i] = rank[i] - mean_y;
            variance_y += diff_y[i]*diff_y[i];
            sum_diff_y += diff_y[i];
        }

        std::vector<std::vector<double>> mean_x = std::vector<std::vector<double>> (m, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> S1 = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));

        #pragma omp parallel
        {
//            std::cout << "computing correlation measure from thread " << omp_get_thread_num() << std::endl;
            std::vector<std::vector<double>> mean_x_local = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));
            std::vector<std::vector<double>> S1_local = std::vector<std::vector<double>>(m, std::vector<double>(n, 0.0));

            int threadnum = omp_get_thread_num();
            int low = sample_size * threadnum / num_thread;
            int high = sample_size * (threadnum+1) / num_thread;
            double ratio = 1.0/sample_size;

            for (int i = low; i < high; ++i){
                for (int j = 0; j < n; ++j){
                    mean_x_local[samples[i][j]][j] += ratio;
                    S1_local[samples[i][j]][j] += diff_y[i];
                }
            }

            #pragma omp critical
            for (int i = 0; i < m; ++i){
                for (int j = 0; j < n; ++j){
                    mean_x[i][j] += mean_x_local[i][j];
                    S1[i][j] += S1_local[i][j];
                }
            }
        }

        std::vector<std::vector<double>> variance_x(m, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> variance_xy(m, std::vector<double>(n, 0.0));
        for (int i = 0; i < m; ++i){
            min_cbm_i[i] = 1.0;
        }
        for (int j = 0; j < n; ++j){
            min_cbm_j[j] = 1.0;
        }
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                variance_x[i][j] = mean_x[i][j]*(1-mean_x[i][j])*sample_size;
                variance_xy[i][j] = (1-mean_x[i][j])*S1[i][j] - mean_x[i][j]*(sum_diff_y - S1[i][j]);
                if (variance_x[i][j] > 0){
                    corr_xr[i][j] = variance_xy[i][j]/sqrt(variance_x[i][j]*variance_y);
                }else if (mean_x[i][j] == 0) {
                    corr_xr[i][j] = 1.0;
                } else{
                    corr_xr[i][j] = -1.0;
                }
                if (corr_xr[i][j] < min_rcm_i[i]){
                    min_rcm_i[i] = corr_xr[i][j];
                }
                if (corr_xr[i][j] < min_rcm_j[j]){
                    min_rcm_j[j] = corr_xr[i][j];
                }
            }
        }
    }

    void Reduce::constructing_test_data(){
        compute_cost_measure();
        compute_correlation_based_measure();
        compute_ranking_based_measure();
        compute_objective_based_measure();
        compute_ranking_correlation_measure();


        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
//            std::cout << "constructing_test_data from thread " << threadnum << std::endl;
            std::string test_s = test_data_dir + g.get_file_name() + test_data_name + std::to_string(threadnum) + ".csv";

            char test_data[test_s.size()+1];
            strcpy(test_data, test_s.c_str());
            std::ofstream test_file(test_data, std::ios::trunc);
            if (! test_file.is_open()){
                std::cout << "Cannot open the output file " <<  test_data << "\n";
            }

            int low = m * threadnum / num_thread;
            int high = m * (threadnum + 1) / num_thread;

            for (int i = low; i < high; ++i){
                for (int j = 0; j < n; ++j){
                    test_file << 0 << " ";
                    test_file << "1:" << std::fixed << std::setprecision(6) << corr_xy[i][j] / min_cbm_i[i] << " ";
                    test_file << "2:" << std::fixed << std::setprecision(6) << corr_xy[i][j] / min_cbm_j[j] << " ";
                    test_file << "3:" << std::fixed << std::setprecision(6) << corr_xr[i][j] / min_rcm_i[i] << " ";
                    test_file << "4:" << std::fixed << std::setprecision(6) << corr_xr[i][j] / min_rcm_j[j] << " ";
                    test_file << "5:" << std::fixed << std::setprecision(6) << ranking_scores[i][j] / max_rbm_i[i] << " ";
                    test_file << "6:" << std::fixed << std::setprecision(6) << ranking_scores[i][j] / max_rbm_j[j] << " ";
                    test_file << "7:" << std::fixed << std::setprecision(6) << objective_scores[i][j] / max_obm_i[i] << " ";
                    test_file << "8:" << std::fixed << std::setprecision(6) << objective_scores[i][j] / max_obm_j[j] << " ";
                    test_file << "9:" << std::fixed << std::setprecision(6) << cost_car[i][j] / max_cost_i[i] << " ";
                    test_file << "10:" << std::fixed << std::setprecision(6) << cost_car[i][j] / max_cost_j[j] << " "<<"\n";
//                    test_file << "7:" << std::fixed << std::setprecision(6) << (double) g.cars_in_class[i] / g.cars << " ";
//                    test_file << "8:" << std::fixed << std::setprecision(6) << g.option_utilisation_class[i]  / g.options  << " ";
//                    test_file << "5:" << std::fixed << std::setprecision(6) << g.max_pq_class[i]  << " ";
//                    test_file << "6:" << std::fixed << std::setprecision(6) << g.min_pq_class[i]  << " ";
//                    test_file << "7:" << std::fixed << std::setprecision(6) << g.ave_pq_class[i]  << " " ;
//                    test_file << "9:" << std::fixed << std::setprecision(6) << double (j+1) / g.cars << " " <<"\n";
                }
            }
            test_file.close();
        }
    }

    void Reduce::loading_output_data(){
        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
//            std::cout << "loading predicted value from thread " << threadnum << std::endl;

            std::string output_s = test_data_dir + g.get_file_name() + output_file_name + std::to_string(threadnum) + ".csv";
            char output_file[output_s.size()+1];
            strcpy(output_file, output_s.c_str());
            std::ifstream predicted_file(output_file);
            if (! predicted_file){
                std::cout << "fail to read the predicted file \n";
            }

            int low = m * threadnum / num_thread;
            int high = m * (threadnum + 1) / num_thread;
            bool value;
            for (int i = low; i < high; ++i){
                for (int j = 0; j < n; ++j){
                    predicted_file >> value;
                    predicted_value[i][j] = value;
                }
            }
            predicted_file.close();
            remove(output_file);
//            const int rem_result = remove(output_file);
//            if(rem_result == 0){
//                std::cout << "Successfully remove predicted data file from thread " << threadnum << std::endl;
//            } else {
//                std::cout << "No such predicted data file from thread " << threadnum << std::endl;
//            }
        }
    }

    // Removing variables using SVM
    void Reduce::removing_variables_svm(){
        std::cout << "problem reduction using SVM" << std::endl;

        constructing_test_data();

        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
//            std::cout << "SVM prediction from thread " << threadnum << std::endl;

            std::string test_s = test_data_dir + g.get_file_name() + test_data_name + std::to_string(threadnum) + ".csv";
            std::string output_s = test_data_dir + g.get_file_name() + output_file_name + std::to_string(threadnum) + ".csv";
            std::string model_s = training_data_dir + training_model_name;
            char test_data[test_s.size()+1];
            char output_file[output_s.size()+1];
            char model_file[model_s.size()+1];
            strcpy(test_data, test_s.c_str());
            strcpy(output_file, output_s.c_str());
            strcpy(model_file, model_s.c_str());

            #pragma omp critical
            svm_predict_model(test_data, model_file, output_file, probability);

//            std::cout << "SVM prediction from thread " << threadnum << " done" << std::endl;
        }

        loading_output_data();
        repair();
        compute_reduced_problem_size();
        std::cout << "problem reduction using SVM done" << std::endl;
    }


    // Removing variables using Linear SVM (fast)
    void Reduce::removing_variables_svm_linear(){
        std::cout << "problem reduction using linear SVM" << std::endl;
        constructing_test_data();

        #pragma omp parallel
        {
            int threadnum = omp_get_thread_num();
//            std::cout << "SVM prediction from thread " << threadnum << std::endl;

            std::string test_s = test_data_dir + g.get_file_name() + test_data_name + std::to_string(threadnum) + ".csv";
            std::string output_s = test_data_dir + g.get_file_name() + output_file_name + std::to_string(threadnum) + ".csv";
            std::string model_s = training_data_dir + training_model_name;
            char test_data[test_s.size()+1];
            char output_file[output_s.size()+1];
            char model_file[model_s.size()+1];
            strcpy(test_data, test_s.c_str());
            strcpy(output_file, output_s.c_str());
            strcpy(model_file, model_s.c_str());

            #pragma omp critical
            linear_svm_predict_model(test_data, model_file, output_file, probability);
//            std::cout << "Linear SVM prediction from thread " << threadnum << " done" << std::endl;

            remove(test_data);
//            const int rem_result = remove(test_data);
//            if(rem_result == 0){
//                std::cout << "Successfully remove test data file from thread " << threadnum << std::endl;
//            } else {
//                std::cout << "No such test data file from thread " << threadnum << std::endl;
//            }
        }

        loading_output_data();
        repair();
        compute_reduced_problem_size();
        std::cout << "Linear SVM prediction done" << std::endl;
    }

    //adding the best solution found in sampling to make sure the graph is still connected
    void Reduce::repair(){
        for (int j = 0; j < n; ++j){
            if (predicted_value[best_sol_found[j]][j] < 0.5){
                predicted_value[best_sol_found[j]][j] = 1;
            }
        }
    }

    void Reduce::compute_reduced_problem_size(){
        int count = 0;
        for (int i = 0; i < m; ++i){
            for (int j = 0; j < n; ++j){
                if(predicted_value[i][j] == 1)
                    count += 1;
            }
        }
        num_variable_left = count;
        std::cout << "proportion of remaining problem size is " << (double)count/(m*n) << std::endl;
    }

    double Reduce::get_wall_clock(){
        struct timeval time;
        if (gettimeofday(&time,NULL)){
            return 0;
        }
        return (double)time.tv_sec + (double)time.tv_usec * .000001;
    }

}
