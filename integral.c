//#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <mpi.h>

#define h 1e-10
#define TO 1.0
#define FROM 0.0

double f (double x)
{
    return x * x * x;
}

double calculate (double l, double r)
{
    long int parts = (r - l) / h;
    double result = f(l) + f(r);
    long int i = 0;

    for ( i = 1; i <= parts; i++ )
    {
        if ( i % 2 == 0 )
            result += 2.0 * f(l + i * h);
        else
            result += 4.0 * f(l + i * h);
    }

    return result / 3.0 * h;
}

int main (int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Status stat;

    int s, r;
    MPI_Comm_size(MPI_COMM_WORLD, &s);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    double segment = (TO - FROM) / s;
    double res = calculate(FROM + r * segment, FROM + (r + 1) * segment);
    //printf("[%lf; %lf]\n", FROM + r * segment, FROM + (r + 1) * segment);

    if ( r == 0 )
    {
        double tmp = 0;

        for(int i = 1; i < s; i++)
        {
            MPI_Recv(&tmp, 1,  MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
            res += tmp;
        }

        printf ("Result is %lf\n", res);
    } else
    {
        MPI_Send(&res, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }


    MPI_Finalize();
    exit(EXIT_SUCCESS);
    return 0;
}
