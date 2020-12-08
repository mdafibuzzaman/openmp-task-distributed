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


	////// create sub matrices in each procs 

	int *sub_mat ; 

	int *sub_vec ; 

	int *total_vec ; 

	int *temp_vec ; 

	int block_size , M , N ;

	int total = atoi(argv[1]) ;  

	M = total ; 
	N = total ; 
	block_size = total / 2 ; 

	int i ,j,k; 

	sub_mat = (int*)malloc( M  * N / block_size * sizeof(int));
	sub_vec = (int*)malloc( block_size * sizeof(int));
	temp_vec = (int*)malloc( block_size * sizeof(int));
	total_vec = (int*)malloc(  N  * sizeof(int));

	#pragma omp parallel
	{
		#pragma omp single 
		{
			#pragma omp task depend (out : sub_mat)
			{
			for (i = 0; i <  M / block_size  ; i++) 
		      for (j = 0; j < N  ; j++) 
		         *(sub_mat + i* N  + j) = rand() % 5 ; 
		 	}

		     #pragma omp task depend (out : sub_vec)
		 		{
		    		for(k = 0 ; k < block_size ; k++){
		    			sub_vec[k] = rand() % 5 ; 
		    			printf("mpi rank %d sub_vec[%d] = %d\n",world_rank,k,sub_vec[k] );
		    		}
		    	}
		   // if(world_rank == 0 )
		    {
		    #pragma omp task depend (in : sub_vec)
		    	{
		    	MPI_Send(sub_vec , block_size , MPI_INT , buddy , 0 , MPI_COMM_WORLD); 
		    	}
		    }
		    //else
		    {
		    #pragma omp task depend (out : temp_vec) 
		    	{
		    		int flag = 0 ; 
		    		MPI_Request req ; 
		    		MPI_Irecv(temp_vec , block_size , MPI_INT , buddy , 0 , MPI_COMM_WORLD , &req);

		    		MPI_Test(&req , &flag , MPI_STATUS_IGNORE) ; 
		    		while (flag == 0) {
		    			#pragma omp taskyield
		    			MPI_Test(&req , &flag , MPI_STATUS_IGNORE) ; 
		    		}

		    		//MPI_Wait(&req , MPI_STATUS_IGNORE);
		    	}
		    }

		    #pragma omp taskwait

		    /*int token;
			if (world_rank != 0) {
    			MPI_Recv(temp_vec, block_size, MPI_INT, world_rank - 1, 0,
             	MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			//printf("Process %d received token %d from process %d\n",
           		//world_rank, token, world_rank - 1);
			} 
			else {
    // Set the token's value if you are process 0
    		token = -1;
			}
			MPI_Send(sub_vec, block_size, MPI_INT, (world_rank + 1) % world_size,
         		0, MPI_COMM_WORLD);

			// Now process 0 can receive from the last process.
			if (world_rank == 0) {
    		MPI_Recv(temp_vec, block_size , MPI_INT, world_size - 1, 0,
             	MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    		//printf("Process %d received token %d from process %d\n",
           	//world_rank, token, world_size - 1);
			}*/

		    #pragma omp task depend(out : total_vec ) depend(in : temp_vec)
			for (i = 0 ; i < block_size ; i++){
				total_vec[buddy*block_size+i] = temp_vec[i];
				//printf("total_vec[%d] = temp_vec[%d] = %d\n",);

			}
		    #pragma omp task depend(out : total_vec ) depend(in : sub_vec)
			for (i = 0 ; i < block_size ; i++){
				total_vec[world_rank*block_size+i] = sub_vec[i];

			}

			///////multiplication////////////








		    #pragma omp taskwait 
		}

    }

    #pragma omp taskwait

    #pragma omp task depend(in : temp_vec)

    for(k = 0 ; k < block_size ; k++){
		    			//temp_vec[k] = rand() % 5 ; 
		    			printf("mpi rank %d temp_vec[%d] = %d\n",world_rank,k,temp_vec[k] );
		    		}


	#pragma omp task depend(in : total_vec)
	for(k = 0 ; k < N ; k++){
		    			//temp_vec[k] = rand() % 5 ; 
		    			printf("mpi rank %d total_vec[%d] = %d\n",world_rank,k,total_vec[k] );
		    		}







	MPI_Finalize();



	return 0;


}