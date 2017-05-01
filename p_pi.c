#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Status status;

    unsigned long int cir, i, max_iter, iter_for_each;
    int s, r;

    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &s);

    double x, y;
    srand(time(NULL));


    max_iter = 4e8;//18e15;

    if ( r == 0 )
        iter_for_each = max_iter / s + max_iter % s;
    else
        iter_for_each = max_iter / s;

    cir = 0;

    for (i=0; i < iter_for_each; i++)
    {
        x = ((double)rand())/RAND_MAX;
        y = ((double)rand())/RAND_MAX;

        if (x*x + y*y <= 1)
        {
            cir++;
        }
    }

    if ( r == 0 )
    {
        unsigned long int res[s];

        for(int j = 1; j < s && s > 1; j++)
        {
            MPI_Recv(res + j,
                     1, MPI_UNSIGNED_LONG,
                     j, 0, MPI_COMM_WORLD, &status);
        }

        unsigned long int sum = cir;

        for(int j = 1; j < s && s > 1; j++)
        {
            sum += res[j];
        }

        long double pi = 4*((long double)sum)/max_iter;
        printf("%Lf\n", pi);

    } else
    {
        MPI_Send(&cir, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
