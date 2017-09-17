#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>


/* Some testing */
int somef(std::vector<int> &data, std::vector<int> &out_data, int node_id);

void ThomasAlgorithm_P(int mynode, int numnodes,
                       int N, double *b, double *a, double *c, double *x, double *q);

void ThomasAlgorithm(int N, double *b, double *a, double *c, double *x, double *q);
void ThomasAlgorithmLU(int N, double *b, double *a, double *c, double *l, double *u, double *d);
void ThomasAlgorithmSolve(int N, double *l, double *u, double *d, double *x, double *q);

//void FillMatrix(int N, double *a, double *b, double *c, double *f, struct bio_params *bio_info);
double CalculateErr(int N, double x);
