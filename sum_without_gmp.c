#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Status status;

    int s, r, n, dn;
    n = 100000;

    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &s);

    if ( r != 0 )
        dn = n / s;
    else
        dn = n / s + n % s;

    double res_f = 0;

    for(int i = 1 + r; i <= n; i+=s)
    {
        double tmp = 1.0;

        for(int j = 2; j <= i; j++)
            tmp = tmp / j;

        res_f += tmp;
    }

    if ( r == 0 )
    {
        for(int j = 1; j < s; j++)
        {
            double tmp_f;
            MPI_Recv(&tmp_f, 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, &status);
            res_f += tmp_f;
        }

        printf("%.19lf\n", res_f);
    } else
    {
        MPI_Send(&res_f, 1, MPI_DOUBLE, 0,  0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
