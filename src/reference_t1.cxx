#include <stdio.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algos.hpp>


#include <BiserLikeModel/biosensor_information.h>
#include <BiserLikeModel/explicit_calculator.h>

using namespace BiserLikeModel;

void callback_crunched(void *ptr, int time) {
    printf("%ds simulated\n", time);
}



void static_fill(struct bio_params *bio_info, int n) {


    // M
    bio_info->km1 = 9.6 * 1e-3;
    bio_info->km2 = 5 * 1e-4;

    // M/s
    bio_info->vmax1 = 1.9 * 1e-4;
    bio_info->vmax2 = 3.9 * 1e-4;


    // [s]
    bio_info->dt = 1e-5;
    bio_info->n = n;
    bio_info->resp_t_meth = FIXED_TIME;

    // [s]
    bio_info->min_t = 100;

    // [s]
    bio_info->resp_t = 20;

    bio_info->dimensionless = false;

    bio_info->out_file_name = "output.dat";
    bio_info->write_to_file = true;

    bio_info->ne = 1;

    // M
    bio_info->pr_0 = 0 * 1e-3;
    bio_info->l_0 = 2e-3;
    bio_info->o2_0 = 2.5 * 1e-4;

    bio_info->rho = 0.56;
    bio_info->layer_count = 3;
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = 1;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dl = 2.2 * 1e-6;
    bio_info->layers[0].Do2 = 0.8 * 1e-5;
    bio_info->layers[0].Dpr = 2.2 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[0].d = 0.025;

    // 1
    bio_info->layers[1].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dl = 6.7 * 1e-6;
    bio_info->layers[1].Do2 = 2.4 * 1e-5;
    bio_info->layers[1].Dpr = 6.7 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[1].d = 0.005;

    // 2
    bio_info->layers[2].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[2].Dl = 6.7 * 1e-6;
    bio_info->layers[2].Do2 = 2.4 * 1e-5;
    bio_info->layers[2].Dpr = 6.7 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[2].d = 0.054;
}


void calc_ref(int n) {
    struct bio_params *bio_info = new bio_params;
    std::vector<double> P, L, t, CP, CL, OP, Chr,  points;

    static_fill(bio_info, n);
    two_layer_model(bio_info, NULL, &callback_crunched, &points, &P, &L, &t, &CL, &CP, &OP, &Chr);

    free(bio_info->layers);
    free(bio_info);
}



int main (int argc, char ** argv) {
    clock_t begin = clock();
    int size_N = 8, N = size_N, i;
    /* Last argument matrix size N */
    if (argc > 1) {
        std::cout << "Parsed parameter N = " << std::stoi(argv[argc-1]) << std::endl;
        size_N = std::stoi(argv[argc-1]);
        N = size_N;
    }

    calc_ref(N);
    clock_t end = clock();

    double time = ((double)(end - begin) / CLOCKS_PER_SEC) * 1000.0;
    std::cout << "sugaista: " << time << std::endl;
    std::ostringstream file_name;
    file_name << "../output/ref_N_"<< size_N << "_p_" << 1;
    std::ofstream file(file_name.str());
    file << time;
    file.close();

    printf("And we done!\n");

    return 0;
}
