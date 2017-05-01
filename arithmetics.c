#include <stdio.h>
#include <mpi.h>
#include <string.h>

#define base 10000
#define number_len 200000

typedef struct
{
    short repr[number_len];
} big_int;

big_int from_string(const char* string)
{
    big_int a = {};
    int len = strlen(string);
    int i = len;
    int p = 0;

    while (i > 0)
    {
        a.repr[p] = 0;

        if ( i >= 4 )
        {
            int mul = 1;
            for(int j = 1; j <= 4; j++)
            {
                a.repr[p] += (string[i-j] - '0') * mul;
                mul *= 10;
            }

            p += 1;
            i -= 4;
        } else
        {
            int mul = 1;
            for(int j = i; j > 0; j--)
            {
                a.repr[p] += (string[j-1] - '0') * mul;
                mul *= 10;
            }

            i = 0;
            p += 1;
        }
    }

    return a;
}

void print(big_int* a)
{
    int flag = 0;

    for( int i = number_len - 1; i >= 0; i-- )
    {
        if ( !flag && a->repr[i] != 0 )
        {
            flag = 1;
            printf("%d", a->repr[i]);
            continue;
        }

        if ( flag )
            printf("%04d", a->repr[i]);
    }

    if ( !flag )
        printf("0");

    printf("\n");
}

big_int sum(int from, int to, big_int* a, big_int* b)
{
    big_int c = {};
    short tmp = 0;
    short r = 0;

    for( int i = from; i < to; i++ )
    {
        tmp = a->repr[i] + b->repr[i];
        c.repr[i] = tmp % base + r;
        r = tmp / base;
    }

    return c;
}

int main(int argc, char** argv)
{
    if ( argc != 3 )
    {
        printf("Enter TWO numbers\n");
        return 0;
    }

    MPI_Init(&argc, &argv);
    MPI_Status status;

    int s, r, from, to;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &s);

    big_int a, b, c0, c1;
    a = from_string(argv[1]);
    b = from_string(argv[2]);

    if ( r == 0 )
    {
        print(&a);
        print(&b);
    }

    from = r * (number_len / s);
    to = (number_len / s) * (r + 1) + (r + 1) / s * (number_len % s);

    if ( r == 0 )
    {
        c0 = sum(from, to, &a, &b); // without carry from previous

        if ( s > 1 )
            MPI_Send(c0.repr + (to - 1), 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

        for( int j = 1; j < s; j++ )
        {
            MPI_Recv(c0.repr + j * (number_len / s),
                     (number_len / s) + (j + 1) / s * (number_len % s),
                     MPI_SHORT,
                     j,
                     0,
                     MPI_COMM_WORLD, &status);
        }

        print(&c0);

    } else
    {
        c0 = sum(from, to, &a, &b); // without carry from previous
        a.repr[from] += 1;
        c1 = sum(from, to, &a, &b); // with carry from previous

        int carry = 0;

        MPI_Recv(&carry, 1, MPI_INT, r - 1, 0, MPI_COMM_WORLD, &status);

        if ( r != s - 1 )
            if ( carry )
                MPI_Send(c1.repr + (to - 1), 1, MPI_INT, r + 1, 0, MPI_COMM_WORLD);
            else
                MPI_Send(c0.repr + (to - 1), 1, MPI_INT, r + 1, 0, MPI_COMM_WORLD);

        if ( carry )
            MPI_Send(c1.repr + from, to - from, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
        else
            MPI_Send(c0.repr + from, to - from, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
