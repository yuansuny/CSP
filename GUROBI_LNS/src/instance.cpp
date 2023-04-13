#include "instance.h"
#include "Random.h"
#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>

using namespace std;

namespace CSP {
    Instance::Instance(std::string file_name, std::string input_dir) :  file_name{file_name}, input_dir{input_dir} {
        cout << file_name << endl;
        read_instance();
        compute_window();
        compute_alpha_beta();
        display();
        read_optimal_solution();
    }

    void Instance::read_instance(){
        string input_file;
        string input_file1 = input_dir  + file_name;

        ifstream file1(input_file1);
        if (file1){
            input_file = input_file1;
        }else{
            input_file = input_file1 + ".txt";
        }

        ifstream file(input_file);
        if (!file){
            std::cout << "problem instance does not exist \n";
            return;
        }

        file >> this->cars;
        file >> this->options;
        file >> this->classes;

        cout << "cars is " << cars << endl;
        cout << "options is " << options << endl;
        cout << "classes is " << classes << endl;

        int t;
        for(int i = 0; i < this->options; i++) {
            file >> t;
            this->max_car_in_option.push_back(t);
        }
        for(int i = 0; i < this->options; i++) {
            file >> t;
            this->max_opt_block_size.push_back(t);
        }
        for(int i = 0; i < this->classes; i++){
            file >> t; 	// id
            file >> t;	// cars in classes
            this->cars_in_class.push_back(t);
            for(int j = 0; j < t; j++) class_of_car.push_back(i);
            car_opt_req.push_back(vector<int>(options,0));
            for(int j = 0; j < this->options; j++){
                file >> t;
                this->car_opt_req[i][j] = t;
                //cout << this->car_opt_req[i][j] << " ";
            }
        }

        //compute pq
        double val;
        for(int i = 0; i < this->options; i++) {
            val = (double) max_car_in_option[i] / max_opt_block_size[i];
            this->pq.push_back(val);
        }

        // compute pq related features
        option_utilisation_class = vector<double>(this->classes, 0.0);
        max_pq_class = vector<double>(this->classes, 0.0);
        min_pq_class = vector<double>(this->classes, 1.0);
        ave_pq_class = vector<double>(this->classes, 0.0);
        for (int i = 0; i < this->classes; ++i){
            for (int j = 0; j < this->options; ++j){
                if (car_opt_req[i][j] > 0.5){
                    option_utilisation_class[i] += 1.0;
                    ave_pq_class[i] += pq[j];
                    if (max_pq_class[i] < pq[j]){
                        max_pq_class[i] = pq[j];
                    }
                    if (min_pq_class[i] > pq[j]){
                        min_pq_class[i] = pq[j];
                    }
                }
            }
        }

        for (int i = 0; i < this->classes; ++i){
            ave_pq_class[i] /= option_utilisation_class[i];
        }
    }

    void Instance::compute_window(){
        int val = 1;
        for(int i = 0; i < this->options; ++i){
            val=lcm(val, this->max_opt_block_size[i]);
        }
        window_size =  val;
    }

    int Instance::lcm(int a, int b){
        int g;
        g=gcd(a,b);
        return a/g*b;
    }

    int Instance::gcd (int a, int b){
        if( a==0 || b==0){
            cout << "Error: option size is 0\n";
            exit(1);
        }
        if (a<b)
            return gcd(b,a);
        else if (a%b== 0)
            return b;
        else
            return gcd(b, a%b);
    }

    void Instance::compute_alpha_beta(){
        Random *r;
        r = new Random(1279948744);
        for(int i = 0; i < cars; i++){
            alpha_weights.push_back(vector<double>(options,0.0));
            beta_weights.push_back(vector<double>(options,0.0));
        }
        for(int i = 0; i < cars; i++){
            for(int j = 0; j < options; j++){
                alpha_weights[i][j] = r->next();
                beta_weights[i][j] = r->next();
            }
        }
        delete r;
    }

    void Instance::display(){
        std::cout << "\nCars: " << this->cars << " ,Options: " << this->options << " ,Classes: " << this->classes << std::endl << "Max cars in options: ";
        for(int i=0;i<this->options;i++) cout << this->max_car_in_option[i] << " ";
        std::cout << "\nMax number block size: " ;
        for(int i=0;i<this->options;i++) std::cout << this->max_opt_block_size[i] << " ";
        std::cout <<  "\nCars in classes: ";
        for(int i=0;i<this->classes;i++) std::cout << this->cars_in_class[i] << " ";
        std::cout <<  "\nClass of car: ";
        for(int i=0;i<this->cars;i++) std::cout << this->class_of_car[i] << " ";
        std::cout <<  "\nOptions required for a class: " << std::endl;
        for(int i=0;i<this->classes;i++){
            for(int j=0;j<this->options;j++) std::cout << this->car_opt_req[i][j] << " ";
            std::cout << std::endl;
        }
        cout << "\nWindow size is "<< this->window_size << "\n";
    }

    void Instance::read_optimal_solution(){
        std::string opt_file_name = input_dir + file_name + ".sol";
        std::ifstream opt_file(opt_file_name);
        if (!opt_file){
            std::cout << "optimal solution is not provided \n";
            return;
        }
        std::string line;
        bool READ_Point = 0;
        int car;
        while(!opt_file.eof()) {
            getline(opt_file, line);
            if (line == "EOF" || line == "-1") {
                break;
            }
            if (READ_Point){
                std::stringstream stream(line);
                while (stream >> car) {
                    opt_solution.push_back(car);
                }
            }
            if (line == "OPTIMAL_SOLUTION"){
                READ_Point = 1;
            }
        }

        std::vector<std::vector<double>> cost_car = std::vector<std::vector<double>>(classes, std::vector<double>(cars, 0.0));
        opt_obj = compute_costs(opt_solution, cost_car);
        cout << "optimal objective value is " << opt_obj << endl;

        // Dubug
        double sum = 0;
        for (int i = 0; i < this->classes; ++i){
            for (int j = 0; j < this->cars; ++j){
                sum += cost_car[i][j];
            }
        }
        cout << "optimal objective value computed by adding each car is " << sum << endl;
//        assert (sum == opt_obj);

        optimal_solution_binary = std::vector<std::vector<bool>>(this->classes, std::vector<bool>(this->cars, 0));
        for (int j = 0; j < this->cars; ++j){
            optimal_solution_binary[opt_solution[j]][j] = 1;
        }

    }

    double Instance::compute_objective_value(const vector<int> &sequence) const {
        double total_uua = 0;
        double total_uoa = 0;
        int l_j_t, min_p, total_options, clas;
        for(int j = 0; j < cars; ++j){
            if(sequence[j] == -1) {
                cout << "Warning, some variables were not bound" << endl;
                return 0;// we are done
            }
            for(int i = 0; i < options; ++i){
                if(j + 1 - max_opt_block_size[i] > 0){
                    l_j_t = j - max_opt_block_size[i] + 1;
                }else{
                    l_j_t = 0;
                }
                if(max_car_in_option[i] < j+1){
                    min_p = max_car_in_option[i];
                }else{
                    min_p = j+1;
                }
                total_options = 0;
                for(int k = l_j_t; k <= j; ++k){ // start from the current car up to the sequence size
                    clas = sequence[k];
                    total_options += car_opt_req[clas][i];
                }
                if(min_p - total_options >= 0){
                    total_uua += double(min_p - total_options) * beta_weights[j][i];
                } else {
                    total_uoa += double(total_options - min_p)* alpha_weights[j][i];
                }
            }
        }
        return total_uua + total_uoa;
    }

    double Instance::compute_costs (const vector<int> &sequence, vector<vector<double>> &cost_car) const {
        double total_uua = 0;
        double total_uoa = 0;
        int l_j_t, min_p, total_options, clas;

//        std::vector<std::vector<double>> cost_car = std::vector<std::vector<double>>(classes, std::vector<double>(cars, 0.0));

        for(int j = 0; j < cars; ++j){
            if(sequence[j] == -1) {
                cout << "Warning, some variables were not bound" << endl;
                return 0;// we are done
            }
            for(int i = 0; i < options; ++i){
                if(j + 1 - max_opt_block_size[i] > 0){
                    l_j_t = j - max_opt_block_size[i] + 1;
                }else{
                    l_j_t = 0;
                }
                if(max_car_in_option[i] < j+1){
                    min_p = max_car_in_option[i];
                }else{
                    min_p = j+1;
                }
                total_options = 0;
                for(int k = l_j_t; k <= j; ++k){ // start from the current car up to the sequence size
                    clas = sequence[k];
                    total_options += car_opt_req[clas][i];
                }
                if(min_p - total_options >= 0){
                    total_uua += double(min_p - total_options) * beta_weights[j][i];
                    // compute the contribution of one car to the objective value
                    for(int k = l_j_t; k <= j; ++k){
                        clas = sequence[k];
                        if (car_opt_req[clas][i] == 0){
                            cost_car[clas][k] += double(min_p - total_options) * beta_weights[j][i] / (j - l_j_t + 1 - total_options);
                        }
                    }
                } else {
                    total_uoa += double(total_options - min_p)* alpha_weights[j][i];
                    // compute the contribution of one car to the objective value
                    for(int k = l_j_t; k <= j; ++k){
                        clas = sequence[k];
                        if (car_opt_req[clas][i] == 1){
                            cost_car[clas][k] += double(total_options - min_p) * alpha_weights[j][i] / total_options;
                        }
                    }

                }
            }
        }
//        cout << "objective value is " << total_uua + total_uoa << endl;
//        double sum = 0;
//        for (int i = 0; i < classes; ++i){
//            for (int j = 0; j < cars; ++j){
//                sum += cost_car[i][j];
//            }
//        }
//        cout << "objective value computed by adding each car is " << sum << endl;

        return total_uua + total_uoa;
    }

//    double Instance::compute_costs (const vector<int> &sequence) const {
//        double total_uua = 0;
//        double total_uoa = 0;
//        int l_j_t, min_p, total_options, clas;
//        for(int i = 0; i < options; ++i){
//            for(int j = 0; j < cars; ++j){
//                if(sequence[j] == -1) {
//                    cout << "Warning, some variables were not bound" << endl;
//                    j = cars; // we are done
//                }
//                else{
//                    if(j - max_opt_block_size[i] > 0){
//                        l_j_t = j - max_opt_block_size[i];
//                    }else{
//                        l_j_t = 0;
//                    }
//                    if(max_car_in_option[i] < j){
//                        min_p = max_car_in_option[i];
//                    }else{
//                        min_p = j;
//                    }
//                    total_options = 0;
//                    for(int k = l_j_t; k < j; ++k){ // start from the current car up to the sequence size
//                        clas = sequence[k];
//                        total_options += car_opt_req[clas][i];
//                    }
//                    if(min_p - total_options > 0){
//                        total_uua += double(min_p - total_options) * beta_weights[j][i];
//                    } else {
//                        total_uoa += double(total_options - min_p)* alpha_weights[j][i];
//                    }
//                }
//            }
//        }
//        return total_uua + total_uoa;
//    }
}