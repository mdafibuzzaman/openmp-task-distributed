#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void trans(int *tot , int *trans , int dim){
	printf("inside trans \n");
	int i , j ; 
	for(j = 0 ; j < dim ; j++){
		for(i = 0 ; i < dim ; i++){
			trans[j*dim+i] = tot[i*dim+j];
		}
	}
}

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

	int M = 4; 
	int N  = 4; 

	int buddy = (world_rank + 1) % world_size ; 

	//printf("rank %d buddy %d\n", world_rank , buddy);


	MPI_Request req ; 


	////// create sub matrices in each procs 

	int *sub_row = (int*)malloc(M * sizeof(int));
	int *total_mat = (int*)malloc(M * N * sizeof(int));
	int *trans_mat = (int*)malloc(M * N * sizeof(int));

	int i ,j,k, flag; 

	
	#pragma omp parallel
	{
		#pragma omp single 
		{
			#pragma omp task depend (out : sub_row) 
			{
			//for (i = 0; i <  M / block_size  ; i++) 
		      for (j = 0; j < M  ; j++) 
		         *(sub_row +  j) = rand() % 5 ; 
		     	//printf("rank %d sub_row[%d] = %d\n", world_rank , j , sub_row[j]);
		 	}


		 	/*#pragma omp task depend(in : sub_row)
		 	{
		 		for (j = 0; j < M  ; j++){
		 			printf("rank %d sub_row[%d] = %d\n", world_rank , j , sub_row[j]);
		 		}	
		 	}*/


		 	#pragma omp task depend(out : total_mat ) depend(in : sub_row)
		 	{
		 		flag = 0 ; 
		 		MPI_Igather(sub_row , M , MPI_INT , total_mat , M , MPI_INT , 2 , MPI_COMM_WORLD , &req);
		 		printf("hello rank %d\n", world_rank);
		 		MPI_Test(&req , &flag , MPI_STATUS_IGNORE) ; 
		    		while (flag == 0) {
		    			#pragma omp taskyield
		    			MPI_Test(&req , &flag , MPI_STATUS_IGNORE) ; 
		    		}


		 	}

		 	//printf("hello rank %d\n", world_rank);

		 	if(world_rank == 2){
		 	#pragma omp task depend(out : trans_mat ) depend(in : total_mat)

		 		{
		 			//trans(total_mat , trans_mat , M);
		 			for(j = 0 ; j < M ; j++){
						for(i = 0 ; i < M ; i++){
							trans_mat[j*M+i] = total_mat[i*M+j];
						}
					}
		 		}


		 	}


			//#pragma omp taskwait

		 	#pragma omp task depend(inout: trans_mat)
		 	{
		 		flag = 0;
		 		MPI_Ibcast(trans_mat , M * M , MPI_INT , 2 , MPI_COMM_WORLD , &req);
		 		MPI_Test(&req , &flag , MPI_STATUS_IGNORE) ; 
		    		while (flag == 0) {
		    			#pragma omp taskyield
		    			MPI_Test(&req , &flag , MPI_STATUS_IGNORE) ; 
		    		}
		 	}






		     
		    	
		    

		    








		    //#pragma omp taskwait 
		}

    }

    //#pragma omp taskwait




	#pragma omp task depend(in : trans_mat)
	for(k = 0 ; k < M * M ; k++){
		    			//temp_vec[k] = rand() % 5 ; 
		    			printf("mpi rank %d trans_mat[%d] = %d\n",world_rank,k,trans_mat[k] );
		    		}







	MPI_Finalize();



	return 0;


}
