#include "instance.h"
#include "reduce.h"
#include "training.h"
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <iomanip>
#include <sstream>

std::string ToString(int value, int digitsCount){
    std::ostringstream os;
    os << setfill('0') << setw(digitsCount) << value;
    return os.str();
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char* argv[]) {
    using namespace CSP;

//    std::vector<std::string> file_name{"4-72", "6-76", "10-93", "16-81", "19-71", "21-90", "26-82", "36-92", "41-66"};
    std::vector<std::string> file_name{"10-93", "16-81", "19-71", "21-90", "26-82", "36-92", "4-72", "41-66", "6-76", "Ranyu10_100_10_25_00", "Ranyu10_100_10_25_01", "Ranyu10_100_10_25_02", "Ranyu10_300_10_25_00", "Ranyu10_300_10_25_01", "Ranyu10_300_10_25_02", "Ranyu10_500_10_25_00", "Ranyu10_500_10_25_01", "Ranyu10_500_10_25_02", "Rlou_100_5_25_00", "Rlou_100_5_25_01", "Rlou_100_5_25_02", "Rlou_300_5_25_00", "Rlou_300_5_25_01", "Rlou_300_5_25_02", "Rlou_500_5_25_00", "Rlou_500_5_25_01", "Rlou_500_5_25_02", "Rnegbhipqvlou10o_100_10_25_00", "Rnegbhipqvlou10o_100_10_25_01", "Rnegbhipqvlou10o_100_10_25_02", "Rnegbhipqvlou10o_300_10_25_00", "Rnegbhipqvlou10o_300_10_25_01", "Rnegbhipqvlou10o_300_10_25_02", "Rnegbhipqvlou10o_500_10_25_00", "Rnegbhipqvlou10o_500_10_25_01", "Rnegbhipqvlou10o_500_10_25_02", "Rnegblopqvlou10o_100_10_25_00", "Rnegblopqvlou10o_100_10_25_01", "Rnegblopqvlou10o_100_10_25_02", "Rnegblopqvlou10o_300_10_25_00", "Rnegblopqvlou10o_300_10_25_01", "Rnegblopqvlou10o_300_10_25_02", "Rnegblopqvlou10o_500_10_25_00", "Rnegblopqvlou10o_500_10_25_01", "Rnegblopqvlou10o_500_10_25_02", "Rnegbvlou10o_100_10_25_00", "Rnegbvlou10o_100_10_25_01", "Rnegbvlou10o_100_10_25_02", "Rnegbvlou10o_300_10_25_00", "Rnegbvlou10o_300_10_25_01", "Rnegbvlou10o_300_10_25_02", "Rnegbvlou10o_500_10_25_00", "Rnegbvlou10o_500_10_25_01", "Rnegbvlou10o_500_10_25_02", "carseq_100_8_20_04", "carseq_100_8_20_06", "carseq_100_8_20_07", "carseq_100_8_20_11", "carseq_100_8_20_13", "carseq_100_8_20_14", "carseq_100_8_20_16", "carseq_100_8_20_18", "carseq_100_8_20_19", "carseq_100_8_20_21", "carseq_100_8_20_22", "carseq_100_8_20_26", "carseq_100_8_20_33", "carseq_100_8_20_35", "carseq_100_8_20_55", "carseq_100_8_20_60", "carseq_100_8_20_62", "carseq_100_8_20_64", "carseq_100_8_20_67", "carseq_100_8_20_71", "carseq_100_8_20_72", "carseq_100_8_20_73", "carseq_100_8_20_77", "carseq_100_8_20_81", "carseq_100_8_20_82", "carseq_100_8_20_85", "carseq_100_8_20_86", "carseq_100_8_20_88", "carseq_100_8_20_90", "carseq_100_8_20_94", "carseq_100_8_20_98", "carseq_100_8_20_99", "carseq_300_8_20_08", "carseq_300_8_20_10", "carseq_300_8_20_14", "carseq_300_8_20_15", "carseq_300_8_20_23", "carseq_300_8_20_25", "carseq_300_8_20_30", "carseq_300_8_20_32", "carseq_300_8_20_36", "carseq_300_8_20_37", "carseq_300_8_20_38", "carseq_300_8_20_40", "carseq_300_8_20_45", "carseq_300_8_20_47", "carseq_300_8_20_53", "carseq_300_8_20_56", "carseq_300_8_20_62", "carseq_300_8_20_76", "carseq_300_8_20_78", "carseq_300_8_20_80", "carseq_300_8_20_81", "carseq_500_8_20_07", "carseq_500_8_20_08", "carseq_500_8_20_09", "carseq_500_8_20_13", "carseq_500_8_20_14", "carseq_500_8_20_15", "carseq_500_8_20_27", "carseq_500_8_20_28", "carseq_500_8_20_30", "carseq_500_8_20_34", "carseq_500_8_20_38", "carseq_500_8_20_39", "carseq_500_8_20_44", "carseq_500_8_20_52", "carseq_500_8_20_53", "carseq_500_8_20_56", "carseq_500_8_20_60", "carseq_500_8_20_62", "carseq_500_8_20_64", "carseq_500_8_20_65", "carseq_500_8_20_66", "carseq_500_8_20_67", "carseq_500_8_20_68", "carseq_500_8_20_72", "carseq_500_8_20_74", "carseq_500_8_20_77", "carseq_500_8_20_79", "carseq_500_8_20_88", "carseq_500_8_20_93", "hipqhiu_100_5_25_00", "hipqhiu_100_5_25_01", "hipqhiu_100_5_25_02", "hipqhiu_300_5_25_00", "hipqhiu_300_5_25_01", "hipqhiu_300_5_25_02", "hipqhiu_500_5_25_00", "hipqhiu_500_5_25_01", "hipqhiu_500_5_25_02", "hipqmedu_100_5_25_00", "hipqmedu_100_5_25_01", "hipqmedu_100_5_25_02", "hipqmedu_300_5_25_00", "hipqmedu_300_5_25_01", "hipqmedu_300_5_25_02", "hipqmedu_500_5_25_00", "hipqmedu_500_5_25_01", "hipqmedu_500_5_25_02", "lopq820_100_8_20_00", "lopq820_100_8_20_01", "lopq820_100_8_20_02", "lopq820_300_8_20_00", "lopq820_300_8_20_01", "lopq820_300_8_20_02", "lopq820_500_8_20_00", "lopq820_500_8_20_01", "lopq820_500_8_20_02", "negbfixedpq_100_5_25_00", "negbfixedpq_100_5_25_01", "negbfixedpq_100_5_25_02", "negbfixedpq_300_5_25_00", "negbfixedpq_300_5_25_01", "negbfixedpq_300_5_25_02", "negbfixedpq_500_5_25_00", "negbfixedpq_500_5_25_01", "negbfixedpq_500_5_25_02", "negbhipqlou_100_5_25_00", "negbhipqlou_100_5_25_01", "negbhipqlou_100_5_25_02", "negbhipqlou_300_5_25_00", "negbhipqlou_300_5_25_01", "negbhipqlou_300_5_25_02", "negbhipqlou_500_5_25_00", "negbhipqlou_500_5_25_01", "negbhipqlou_500_5_25_02", "negbhiu_100_5_25_00", "negbhiu_100_5_25_01", "negbhiu_100_5_25_02", "negbhiu_300_5_25_00", "negbhiu_300_5_25_01", "negbhiu_300_5_25_02", "negbhiu_500_5_25_00", "negbhiu_500_5_25_01", "negbhiu_500_5_25_02", "nobhiu_100_5_25_00", "nobhiu_100_5_25_01", "nobhiu_100_5_25_02", "nobhiu_300_5_25_00", "nobhiu_300_5_25_01", "nobhiu_300_5_25_02", "nobhiu_500_5_25_00", "nobhiu_500_5_25_01", "nobhiu_500_5_25_02", "pb_200_01", "pb_200_02", "pb_200_03", "pb_200_04", "pb_200_05", "pb_200_06", "pb_200_07", "pb_200_08", "pb_200_09", "pb_200_10", "pb_300_01", "pb_300_02", "pb_300_03", "pb_300_04", "pb_300_05", "pb_300_06", "pb_300_07", "pb_300_08", "pb_300_09", "pb_300_10", "pb_400_01", "pb_400_02", "pb_400_03", "pb_400_04", "pb_400_05", "pb_400_06", "pb_400_07", "pb_400_08", "pb_400_09", "pb_400_10", "posbhiu_100_5_25_00", "posbhiu_100_5_25_01", "posbhiu_100_5_25_02", "posbhiu_300_5_25_00", "posbhiu_300_5_25_01", "posbhiu_300_5_25_02", "posbhiu_500_5_25_00", "posbhiu_500_5_25_01", "posbhiu_500_5_25_02", "randN_100_5_25_00", "randN_100_5_25_01", "randN_100_5_25_02", "randN_300_5_25_00", "randN_300_5_25_01", "randN_300_5_25_02", "randN_500_5_25_00", "randN_500_5_25_01", "randN_500_5_25_02"};

//    std::vector<std::string> file_name{"pb_200_01", "pb_200_02", "pb_200_03", "pb_200_04", "pb_200_05", "pb_200_06", \
//        "pb_200_07", "pb_200_08", "pb_200_09", "pb_200_10", "pb_300_01", "pb_300_02", "pb_300_03", "pb_300_04", \
//        "pb_300_05", "pb_300_06", "pb_300_07", "pb_300_08", "pb_300_09", "pb_300_10", "pb_400_01", "pb_400_02", \
//        "pb_400_03", "pb_400_04", "pb_400_05", "pb_400_06", "pb_400_07", "pb_400_08", "pb_400_09", "pb_400_10"};

    std::vector<std::string> training_file_name{"4-72", "6-76", "21-90", "26-82", "36-92", "41-66"};

//    const std::string input_dir = "../../cs_inst/";
    const std::string input_dir = "../../testset/";
//    const std::string input_dir = "../../cs_inst/ProblemDataSet200to400/";
    const std::string output_dir = "../results/";

    double time_limit = 3600;

    std::uint32_t runs = 1u;
    std::uint32_t d;
    sscanf(argv[1], "%d", &d);

    std::string input_file_name = file_name[d];


    //uncomment this if training is required
//    auto training = Training(training_file_name, "../../cs_inst/");
//    training.generate_training_model_svm();
//    training.generate_training_model_svm_linear();

    std::string output_obj_filename, output_sol_filename, output_time_filename;
    output_obj_filename = output_dir + input_file_name + "_res_obj.csv";
//    output_sol_filename = output_dir + input_file_name + "_res_sol.csv";
    output_time_filename = output_dir + input_file_name + "_res_time.csv";
//    std::ofstream output_file_obj (output_obj_filename, std::ios::trunc);
//    std::ofstream output_file_time (output_time_filename, std::ios::trunc);
    std::ofstream output_file_obj (output_obj_filename, std::ios::app);
    std::ofstream output_file_time (output_time_filename, std::ios::app);
//    std::ofstream output_file_sol (output_sol_filename, std::ios::trunc);

   if (output_file_obj.is_open()){
//        output_file_obj << "Y" << "\n";
        output_file_obj.close();
    } else{
        std::cout << "Cannot open the output file " + output_obj_filename << "\n";
        return 0;
    }
//    output_file_obj.close();

    if (output_file_time.is_open()){
//        output_file_time << "t_total" << "\n";
        output_file_time.close();
    } else{
        std::cout << "Cannot open the output file " + output_time_filename << "\n";
        return 0;
    }
//    output_file_time.close();

//    if (output_file_sol.is_open()){
//        output_file_sol << "OPTIMAL_SOLUTION" << "\n";
//    } else{
//        std::cout << "Cannot open the output file " + output_sol_filename << "\n";
//        return 0;
//    }
//    output_file_sol.close();

    for (auto i = 0u; i < runs; ++i){

        double w0 = get_wall_time();
        double c0 = get_cpu_time();
        const auto instance = Instance(input_file_name, input_dir);

        auto reduce = Reduce(instance);
        double w1 = get_wall_time();
//        double c1 = get_cpu_time();

        reduce.solver_LNS_GUROBI(time_limit - w1 + w0);
//        reduce.solver_LNS_sampling();
        double w2 = get_wall_time();
        double c2 = get_cpu_time();
        std::cout << "Best objective value found is  " << reduce.get_objective_value_found() << "\n";
        std::cout << "Solving wall time is  " << w2 - w0 << "s\n";
        std::cout << "Solving CPU time is  " << c2 - c0 << "s\n";

        output_file_obj.open(output_obj_filename, std::ios::app);
        output_file_obj << reduce.get_objective_value_found() <<"\n";
        output_file_obj.close();

        output_file_time.open(output_time_filename, std::ios::app);
        output_file_time << w2 - w0 << "\n";
        output_file_time.close();
    }

    return 0;
}
