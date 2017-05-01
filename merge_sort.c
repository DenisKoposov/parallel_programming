#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define MAX 100

int* generate_data(int n);
int* merge(int* a1, int n1, int* a2, int n2);
int* merge_sort(int* a1, int l, int );
int* parallel_merge(int* data, int size, \
                    int height, int rank, int num_proc);
int check_answer(int* data, int* check, int n);

int* generate_data(int n)
{
    srand(time(NULL));
    int* data = (int*) malloc (n * sizeof(int));
    
    for(int i = 0; i < n; i++)
        data[i] = rand() % MAX;
        
    return data;
}

int* merge(int* a1, int n1, int* a2, int n2)
{
    int* res = (int*) malloc (sizeof(int) * (n1 + n2));
    int i = 0, j = 0;
    
    while (i < n1 && j < n2)
        res[i+j] = a1[i] <= a2[j] ? a1[i++] : a2[j++];

    while ( i < n1 )
    {
        res[i+j] = a1[i];
        ++i;
    }
    
    while ( j < n2 )
    {
        res[i+j] = a2[j];
        ++j;
    }
    
    return res;
}

int* merge_sort(int* a1, int l, int r)
{
    int n = r - l + 1;

    if ( n <= 1 )
        return a1 + l; 

    int mid = (r + l) / 2;

    int* sorted_left  = merge_sort(a1, l, mid);
    int* sorted_right = merge_sort(a1, mid+1, r);

    int* ans = merge(sorted_left, mid-l+1, \
                     sorted_right, r-mid);
    
    return ans;
}

void parallelMerge(int* data, int size, int height, int rank, int num_proc)
{
    int parent = rank & ~(1 << height);
    int next = height - 1;
    int child = rank | (1 << next);

    //printf("%d %d %d %d\n", parent, height, next, child);

    if (height > 0)
    {
        if (child >= num_proc)
        {
            parallelMerge(data, size, next, rank, num_proc);
        } else
        {
            int left_size = size / 2;
            int right_size = size - left_size;

            int* left_array = (int*) malloc (left_size * sizeof(int));	
            int* right_array = (int*) malloc (right_size * sizeof(int));

            memcpy(left_array, data, left_size * sizeof(int));
            memcpy(right_array, data + left_size, right_size * sizeof(int));

            int info[2];
            info[0] = right_size;
            info[1] = next;
            // sending right part to the child
            MPI_Send(info, 2, MPI_INT, child, 0, MPI_COMM_WORLD);
            MPI_Send(right_array, right_size, MPI_INT, child, 0, MPI_COMM_WORLD);
            // applying algorithm to the left part
            parallelMerge(left_array, left_size, next, rank, num_proc);
            // getting sorted right part from the child
            MPI_Recv(right_array, right_size, MPI_INT, child, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // merging two parts
            int* tmp = merge(left_array, left_size, right_array, right_size);
            memcpy(data, tmp, size);
            
            free(tmp);
            tmp = NULL;

            free(left_array);
            left_array = NULL;

            free(right_array);
            right_array = NULL;
        }
    } else
    {
        int *tmp = (int*) malloc (size * sizeof(int));
        tmp = merge_sort(data, 0, size-1);
        memcpy(data, tmp, size * sizeof(int));
    }

    if (parent != rank)
    {
        MPI_Send(data, size,  MPI_INT, parent, 0, MPI_COMM_WORLD);
    }
}

int check_answer(int* data, int* example, int n)
{
    return !!(memcmp(data, example, n));
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int N = strtol(argv[1], NULL, 10);
    int s, r;

    MPI_Comm_size(MPI_COMM_WORLD, &s);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Barrier(MPI_COMM_WORLD); 
    double start_time;
    
    if (r == 0)
    {
        start_time = MPI_Wtime();

        int* data = generate_data(N);
        int* check = (int*) malloc (N * sizeof(int));
        memcpy(check, data, N * sizeof(int));

        int root_Ht = (int) ceil(log2(s));
     //   printf("root height=%d\n", root_Ht);

        parallelMerge(data, N, root_Ht, r, s);

        if (r == 0)
            printf("Time: %lf\n",  MPI_Wtime() - start_time);

        if (check_answer(data, check, N))
            printf("OK\n");
        else
            printf("FAILED\n");

        free(data);
        data = NULL;
        free(check);
        data = NULL;

        MPI_Finalize();
    } else
    {
        int info[2];
        MPI_Recv(info, 2, MPI_INT, MPI_ANY_SOURCE, \
                 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int *received_data = (int*)calloc(info[0], sizeof(int));
        MPI_Recv(received_data, info[0], MPI_INT, MPI_ANY_SOURCE, \
                 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        parallelMerge(received_data, info[0], info[1], r, s);

        free(received_data);
        received_data = NULL;
//        printf("%d finished\n", r);
        MPI_Finalize();
    }

    return 0;
}
