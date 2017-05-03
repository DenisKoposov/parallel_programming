#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

/*
    Define left part of your equation here:
    u'_t + a * u'_x = f(x, y, y') = t + x
    u(0, x) = x, 0 <= x <= 1
    u(t, 0) = t, 0 <= t <= 1
*/
static inline double f(double x, double t){ return t + x; }
/*
    u_l_k = u(t0 + l*tau, x0 + k*h)
    u_l_l_1 = u(t0 + l*tau, x0 + (k - 1)*h )
    u_next = u(t0 + (l + 1)*tau, x0 + k*h)
*/
static inline double u_next(double x, double t, double u_l_k, \
              double u_l_k_1, double tau, double h, double f, double a)
{
    return u_l_k + tau * f - a * tau / h * (u_l_k - u_l_k_1);
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        printf("Incorrect number of arguments\n");
        exit(EXIT_FAILURE);
    }
    
    MPI_Init(&argc, &argv);
    int K = strtol(argv[1], NULL, 10);
    int L = strtol(argv[2], NULL, 10);
    int rank, size;
    double start_time;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double t0 = 0;
    double tL = 1.0;
    double x0 = 0;
    double xK = 1.0;
    double a = 10;
    double tau = (tL - t0) / L;
    double h = (xK - x0) / K;
    double** u;
    double shift = 0;
    int width = K / size;
    int width_root = width + K % size + 1;
    int width_common = width;

    if ( rank == 0 )
    {
        width = width_root;
        printf("rank:%d K:%d h:%lf L:%d tau:%lf width:%d\n", rank, K, h, L, tau, width);
    } else
    {
        shift = h * width_root + (rank - 1) * h * width;
        printf("rank:%d K:%d h:%lf L:%d tau:%lf width:%d\n", rank, K, h, L, tau, width);
    }

    if ( rank == 0 )
    {
        start_time = MPI_Wtime();
        u = (double**) malloc (sizeof(double*) * 2);
        
        for(int i = 0; i < 2; ++i)
            u[i] = (double*) malloc (sizeof(double) * (K + 1));

    } else
    {
        u = (double**) malloc (sizeof(double*) * 2);
        
        for(int i = 0; i < 2; ++i)
            u[i] = (double*) malloc (sizeof(double) * width);
    }

    // Initial conditions
    for(int i = 0; i < width; ++i)
        u[0][i] = i * h + shift;

    for(int i = 0; i < L; ++i)
    {
        for(int j = 1; j < width; j++)
            u[1][j] = u_next(h * j + shift, tau * i, u[0][j], u[0][j-1], \
                               tau, h, f(h * j + shift, tau * i), a);
        if ( rank == 0 )
        {
            u[1][0] = (i + 1) * tau;
            
            if ( size != 1 )
            {
                double tmp = u[0][width-1]; 
                MPI_Send(&tmp, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            }
        }

        if ( size != 1 && rank != 0 && rank != size - 1 )
        {
            double tmp_recv = 0;
            MPI_Recv(&tmp_recv, 1, MPI_DOUBLE, rank - 1, 0, \
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            double tmp_send = u[0][width-1];
            MPI_Send(&tmp_send, 1, MPI_DOUBLE, \
                     rank + 1, 0, MPI_COMM_WORLD);
        
            u[1][0] = u_next(shift, tau * i, u[0][0], tmp_recv, \
                               tau, h, f(shift, tau * i), a);
        }

        if ( size != 1 && rank == size - 1 )
        {
            double tmp_recv = 0;
            MPI_Recv(&tmp_recv, 1, MPI_DOUBLE, rank - 1, 0, \
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            u[1][0] = u_next(shift, tau * i, u[0][0], tmp_recv, \
                               tau, h, f(shift, tau * i), a);
        }

        double* tmp_pointer = u[0];
        u[0] = u[1];
        u[1] = tmp_pointer;
    }

    // End of computational part. Collecting information
    if ( rank == 0 )
    {
        printf("Time: %lf\n", MPI_Wtime() - start_time);

        for(int i = 1; i < size; ++i)
        {
            MPI_Recv(u[0] + width + (i - 1) * width_common, \
                     width_common, MPI_DOUBLE, i, 0, \
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
        FILE* out = fopen("output.csv", "w");
        
        for(int i = 0; i <= K; ++i)
        {
            fprintf(out, "%lf\n", u[0][i]);
        }
    }

    if ( rank != 0 )
    {
        for(int i = 0; i < width; ++i)
            MPI_Send(u[0], \
                     width, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
