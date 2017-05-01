#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gmp.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Status status;

    int s, r, n, dn;
    n = 40000;

    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &s);

    if ( r != 0 )
        dn = n / s;
    else
        dn = n / s + n % s;

    mpz_t fact;
    mpf_t res_f;
    mpf_t tmp_f;

    mpf_set_default_prec(256);

    mpf_init(res_f);
    mpf_init(tmp_f);
    mpz_init(fact);

    mpf_set_ui(res_f, 0);

    for (unsigned int i = 1 + r; i <=n; i += s)
    {
        mpz_fac_ui(fact, i);
        mpf_set_z(tmp_f, fact);
        mpf_ui_div(tmp_f, 1, tmp_f);
        mpf_add(res_f, res_f, tmp_f);
    }

    if ( r == 0 )
    {
        for(int j = 1; j < s; j++)
        {
            MPI_Recv(tmp_f, sizeof(tmp_f), MPI_BYTE, j, 0, MPI_COMM_WORLD, &status);
            mpf_add(res_f, res_f, tmp_f);
        }

        mpf_out_str(stdout, 10, 0, res_f);
        printf("\n");
    } else
    {
        MPI_Send(res_f, sizeof(res_f), MPI_BYTE, 0,  0, MPI_COMM_WORLD);
    }

    mpz_clear(fact);
    mpf_clear(res_f);
    mpf_clear(tmp_f);

    MPI_Finalize();

    return 0;
}
