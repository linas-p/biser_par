#include <mpi.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <unistd.h>
int fake_time = 30000; // like 3 sec

#include <algos.hpp>

#include <BiserLikeModel/biosensor_information.h>
#include <BiserLikeModel/explicit_calculator.h>

using namespace BiserLikeModel;

void static_fill(struct bio_params *bio_info) {


    // M
    bio_info->km1 = 9.6 * 1e-3;
    bio_info->km2 = 5 * 1e-4;

    // M/s
    bio_info->vmax1 = 1.9 * 1e-4;
    bio_info->vmax2 = 3.9 * 1e-4;


    // [s]
    bio_info->dt = 1e-2;
    bio_info->n = 4;
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
    bio_info->l_0 = 1.0;
    bio_info->o2_0 = 1.0;

    bio_info->rho = 0.5;
    bio_info->layer_count = 3;
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = 1;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dl = 1.0;
    bio_info->layers[0].Do2 = 1.0;
    bio_info->layers[0].Dpr = 1.0;
    // [um] -> [cm]
    bio_info->layers[0].d = 1;

    // 1
    bio_info->layers[1].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dl = 1.0;
    bio_info->layers[1].Do2 = 1.0;
    bio_info->layers[1].Dpr = 1.0;
    // [um] -> [cm]
    bio_info->layers[1].d = 1;

    // 2
    bio_info->layers[2].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[2].Dl = 1.0;
    bio_info->layers[2].Do2 = 1.0;
    bio_info->layers[2].Dpr = 1.0;
    // [um] -> [cm]
    bio_info->layers[2].d = 1;
}


void FillMatrix(int N, double **a, double **b, double **c, double **f,
        double *x_old, struct bio_params *bio_info) {
    *a = new double[N];
    *b = new double[N];
    *c = new double[N];
    *f = new double[N];
    int N_0 = 0, N_Rm = (N/2)-1, N_Rp = (N/2), N_1 = N-1;
    int n = (N-2)/2;
    int i, j;




    double *idx = new double[N];
    for(i = 0; i < N; i++) { idx[i] = i;}

    double dr = bio_info->layers[0].d/n;
    double dr2 = pow(dr, 2);

    for (i = 0; i < N_Rp; i++) {
        (*a)[i] = -2.0/dr2 - 1.0/(1+x_old[i]);
        (*b)[i+1] = 1.0/dr2 - 1/dr * 1/(idx[i+1]*dr);
        (*c)[i+1] = 1.0/dr2 + 1/dr * 1/(idx[i+1]*dr);
        (*f)[i] = x_old[i];
    }
    (*c)[0] = -(*a)[0];
    (*f)[0] = (*f)[0]/2.0;
    (*f)[N_Rm] += ((*b)[N_Rm] + (*a)[N_Rm])*bio_info->l_0;

    double dr_2 = bio_info->layers[1].d/n;
    double dr_2_2 = pow(dr_2, 2);
    for (i = N_Rp; i < N; i++) {
        (*a)[i] = -2.0/dr_2_2;
        (*b)[i+1] = 1.0/dr_2_2 - 1/dr_2 * 1/(idx[i+1]*dr_2);
        (*c)[i+1] = 1.0/dr_2_2 + 1/dr_2 * 1/(idx[i+1]*dr_2);
        (*f)[i] = x_old[i];
    }

    (*f)[N_1] += bio_info->l_0*((*a)[N_1] + (*b)[N_1]);

    /*for(i = 0; i < N; i++) {
        printf("%d %lf %lf %lf %lf\n", i, (*b)[i], (*a)[i],  (*c)[i], (*f)[i]);
    }*/

    delete [] idx;

}


int main(int argc, char ** argv) {

    double starttime, endtime;
    int size_N = 8, N = size_N, i;
    /* Last argument matrix size N */
    if (argc > 1) {
        std::cout << "Parsed parameter N=" << std::stoi(argv[argc-1]) << std::endl;
        size_N = std::stoi(argv[argc-1]);
        N = size_N;
    }

    // std::cout << "Given n = " << size_N << ", " << std::endl;
    struct bio_params *bio_info = new bio_params;
    static_fill(bio_info);

    int mynode, totalnodes;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

    double *a,*b,*c,*q,*x;
    x = new double[N];
    for(i=0; i < N/2; i++) {
        x[i] = 0.5;
    }
    for(i=N/2; i < N; i++) {
        x[i] = 1.0;
    }

    FillMatrix(N, &a, &b, &c, &q, x, bio_info);

    starttime = MPI_Wtime();
    if (totalnodes != 1) {
        ThomasAlgorithm_P(mynode, totalnodes, N,b,a,c,x,q);
    } else {
        ThomasAlgorithm(N,b,a,c,x,q);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    endtime = MPI_Wtime();

    if (mynode == 0) {
        /*for(i=0; i<N; i++)
            std::cout << x[0] << std::endl;
        std::cout << "Done! " << N << std::endl;*/
    }

    if (mynode == 0) {
        double time = (endtime - starttime) * 1000.0;
        std::cout << "sugaista: " << time << std::endl;
        std::ostringstream file_name;
        file_name << "../output/N_"<< size_N << "_p_" << totalnodes;
        std::ofstream file(file_name.str());
        file << time;
        file.close();
    }

    free(bio_info->layers);
    free(bio_info);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] q;
    delete[] x;

    return 0;
}
