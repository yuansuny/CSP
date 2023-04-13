#include "training.h"

namespace CSP {
    Training::Training(std::vector<std::string> training_files, std::string input_dir) : training_files{training_files}, input_dir{input_dir}{
        std::cout << "number of training graphs is: " << training_files.size() << "\n";
        construct_training_set();
    }

    void Training::construct_training_set(){
        std::string train_s = train_data_dir  + train_file_name;
        char train_data[train_s.size()+1];
        strcpy(train_data, train_s.c_str());

        std::ofstream train_file(train_data, std::ios::trunc);
        if (! train_file.is_open()){
            std::cout << "Cannot open the output file " <<  train_data << "\n";
            return;
        }
        train_file.close();

        std::uint64_t num0 = 0;
        std::uint64_t num1 = 0;
        for (int d = 0; d < training_files.size(); ++d){
            const auto ins = Instance(training_files[d], input_dir);
            auto reduce = Reduce(ins);
            reduce.multi_thread_random_sampling();
            reduce.compute_cost_measure();
            reduce.compute_correlation_based_measure();
            reduce.compute_ranking_based_measure();
            reduce.compute_objective_based_measure();
            reduce.compute_ranking_correlation_measure();
            train_file.open(train_data, std::ios::app);
            for (int i = 0; i < ins.classes; ++i){
                for (int j = 0; j < ins.cars; ++j){
                    if (ins.get_optimal_value(i,j) == 1)
                        num1++;
                    else
                        num0++;
                    train_file << ins.get_optimal_value(i,j) << " ";
                    train_file << "1:" << std::fixed << std::setprecision(6) << reduce.get_cbm_value(i,j) / reduce.get_min_cbm_i(i) << " ";
                    train_file << "2:" << std::fixed << std::setprecision(6) << reduce.get_cbm_value(i,j) / reduce.get_min_cbm_j(j) << " ";
                    train_file << "3:" << std::fixed << std::setprecision(6) << reduce.get_rcm_value(i,j) / reduce.get_min_rcm_i(i) << " ";
                    train_file << "4:" << std::fixed << std::setprecision(6) << reduce.get_rcm_value(i,j) / reduce.get_min_rcm_j(j) << " ";
                    train_file << "5:" << std::fixed << std::setprecision(6) << reduce.get_rbm_value(i,j) / reduce.get_max_rbm_i(i) << " ";
                    train_file << "6:" << std::fixed << std::setprecision(6) << reduce.get_rbm_value(i,j) / reduce.get_max_rbm_j(j) << " ";
                    train_file << "7:" << std::fixed << std::setprecision(6) << reduce.get_obm_value(i,j) / reduce.get_max_obm_i(i) << " ";
                    train_file << "8:" << std::fixed << std::setprecision(6) << reduce.get_obm_value(i,j) / reduce.get_max_obm_j(j) << " ";
                    train_file << "9:" << std::fixed << std::setprecision(6) << reduce.get_cost_car(i,j) / reduce.get_max_cost_i(i) << " ";
                    train_file << "10:" << std::fixed << std::setprecision(6) << reduce.get_cost_car(i,j) / reduce.get_max_cost_j(j) << " "<<"\n";
//                    train_file << "7:" << std::fixed << std::setprecision(6) << (double) ins.cars_in_class[i] / ins.cars << " ";
//                    train_file << "8:" << std::fixed << std::setprecision(6) << ins.option_utilisation_class[i]  / ins.options  << " ";
//                    train_file << "5:" << std::fixed << std::setprecision(6) << ins.max_pq_class[i]  << " ";
//                    train_file << "6:" << std::fixed << std::setprecision(6) << ins.min_pq_class[i]  << " ";
//                    train_file << "7:" << std::fixed << std::setprecision(6) << ins.ave_pq_class[i]  << " ";
//                    train_file << "9:" << std::fixed << std::setprecision(6) << double (j+1) / ins.cars << " " <<"\n";
                }
            }
            train_file.close();
        }
        std::cout << "num0 is " << num0 << "; " << "num1 is " << num1 <<  std::endl;
        weight = (double)alpha * num0/num1;
    }

    void Training::generate_training_model_svm(){
        std::string train_s = train_data_dir  + train_file_name;
        std::string model_s = train_data_dir + train_model_name;
        char train_data[train_s.size()+1];
        char model_file[model_s.size()+1];
        strcpy(train_data, train_s.c_str());
        strcpy(model_file, model_s.c_str());
        std::cout << train_data << std::endl;
        std::cout << model_file << std::endl;
        std::cout << weight << std::endl;
        std::cout << kernel_type << std::endl;
        std::cout << prob << std::endl;
        svm_train_model(train_data, model_file, weight, kernel_type, prob);

        const int rem_result = remove(train_data);
        if(rem_result == 0){
            std::cout << "Successfully remove training data file" << std::endl;
        } else {
            std::cout << "No such training data file " << std::endl;
        }
    }

    void Training::generate_training_model_svm_linear(){
        std::string train_s = train_data_dir  + train_file_name;
        std::string model_s = train_data_dir + train_model_name;
        char train_data[train_s.size()+1];
        char model_file[model_s.size()+1];
        strcpy(train_data, train_s.c_str());
        strcpy(model_file, model_s.c_str());
        linear_svm_train_model(train_data, model_file, weight);

        const int rem_result = remove(train_data);
        if(rem_result == 0){
            std::cout << "Successfully remove training data file" << std::endl;
        } else {
            std::cout << "No such training data file " << std::endl;
        }
    }
}
