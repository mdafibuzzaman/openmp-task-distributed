#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc , char	*argv[]){

	int provided;

	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
    if(provided < MPI_THREAD_MULTIPLE)
    {
        printf("The threading support level is lesser than that demanded.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    else
    {
        printf("The threading support level corresponds to that demanded.\n");
    }
	int world_size ;

	MPI_Comm_size(MPI_COMM_WORLD,&world_size);

	int world_rank ; 
	MPI_Comm_rank(MPI_COMM_WORLD , &world_rank ) ;

	//int num_threads = atoi(argv[1]);
	int num_threads = 2 ; 
	omp_set_num_threads(num_threads) ; 

	int omp_thread_cnt , omp_thread_num ; 

	/*#pragma omp parallel default(shared) private(omp_thread_num, omp_thread_cnt)
	{
		omp_thread_cnt = omp_get_num_threads();
		omp_thread_num = omp_get_thread_num();

		printf("omp thread %d out of %d threads mpi rank %d out of %d procs\n",omp_thread_num,omp_thread_cnt, world_rank, world_size);

	} */

	int a ; 
	int b ; 

	int buddy = (world_rank + 1) % world_size ; 

	printf("rank %d buddy %d\n", world_rank , buddy);

#pragma omp parallel num_threads(num_threads)

	{

	#pragma omp single nowait
	{

	if(world_rank == 0)
	{
		#pragma omp task depend(out:a) priority(1)
		a = 10 ; 
 
		#pragma omp task depend(in: a) priority(1)

		MPI_Send(&a , 1 , MPI_INT , buddy , 0 , MPI_COMM_WORLD); 

		#pragma omp task depend(out:a) priority(1)
		a = 20 ; 

	}

	else
	{
		#pragma omp task depend(out:b) priority(2)

		MPI_Recv(&a , 1 , MPI_INT , buddy , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE) ; 
	}
}

}

	#pragma omp taskwait

	/*#pragma omp parallel default(shared) private(omp_thread_num, omp_thread_cnt)
	{
		omp_thread_cnt = omp_get_num_threads();
		omp_thread_num = omp_get_thread_num();

		printf("a = %d omp thread %d out of %d threads mpi rank %d out of %d procs\n",a, omp_thread_num,omp_thread_cnt, world_rank, world_size);

	} */

	printf("a = %d omp thread %d out of %d threads mpi rank %d out of %d procs\n",a, omp_thread_num,omp_thread_cnt, world_rank, world_size);


	MPI_Finalize();



	return 0;


}