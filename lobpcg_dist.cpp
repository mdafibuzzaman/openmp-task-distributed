#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <string.h>

using namespace std ;  





struct block
{
    int nnz;
    int roffset, coffset;
    unsigned short int *rloc, *cloc;
    //int *rloc , *cloc ; 
    double *val;
};

int numcols, numrows, nnonzero, nthrds;
int nrows, ncols, nnz;
int wblk, nrowblks, ncolblks;

int *colptrs, *irem;
double *xrem;


int *nnzPerRow;
block *matrixBlock;



/////////mfdn variables/////////////

int ndiag , nsegments , npe; 
int world_rank , world_size ; 

MPI_Comm col_comm , row_comm , diag_comm ; 

int mysegment ; 


int lvecdim , totaldim ; 


/////////////////mfdn memory variables/////////
int *ilvecshare , *ilvecdisp ; 
double *distamp , *distHamp ; 
double * amp , *Hammp , *ampT , *HammpT ;

int idiag ; 

int ncol , nrow ; 







void read_custom(char* filename, double *&xrem)
{


    FILE *fp;

    double tstart, tend;

    tstart = omp_get_wtime();

    fp = fopen(filename, "rb");

        if(fp == NULL)
    {
            cout << "invalid matrix file name" << endl;
        return;
    }

    fread(&numrows, sizeof(int), 1, fp);
    cout<<"row: "<<numrows<<endl;

    fread(&numcols, sizeof(int), 1, fp);
    cout<<"colum: "<<numcols<<endl;

    fread(&nnonzero, sizeof(int), 1, fp);
    cout<<"non zero: "<<nnonzero<<endl;

    colptrs = new int[numcols + 1];
        irem = new int[nnonzero];
        //xrem = new T[nnonzero];
        xrem = new double[nnonzero];
        float *txrem = new float[nnonzero];
    cout << "Memory allocation finished" << endl;

    fread(colptrs, sizeof(int), numcols+1, fp);
        cout << "finished reading colptrs" << endl;

    fread(irem, sizeof(int), nnonzero, fp);
    cout << "finished reading irem" << endl;

    fread(txrem, sizeof(float), nnonzero, fp);
        cout << "finished reading xrem" << endl;
    for(int i = 0 ; i < nnonzero ; i++){
    
        xrem[i] = txrem[i];
    }


    for(int i = 0 ; i < numcols+1 ; i++)
    {	
    	colptrs[i]--;
    }

    for(int i = 0 ; i < nnonzero ; i++)
    	irem[i]--;

    //for (int i = 0; i < nnonzero; ++i)
    //{
    	
   // 	xrem[i]--;
    //}
    
    delete []txrem;
    tend = omp_get_wtime();
    cout << "Matrix is read in " << tend - tstart << " seconds." << endl;


	

}




void csc2blkcoord(block *&matrixBlock, double *xrem)
{
    int i, j, r, c, k, k1, k2, blkr, blkc, tmp;
    int **top;
    nrowblks = ceil(numrows / (double)(wblk));
    ncolblks = ceil(numcols / (double)(wblk));
    cout << "wblk = " << wblk << endl;
    cout << "nrowblks = " << nrowblks << endl;
    cout << "ncolblks = " << ncolblks << endl;

    //matrixBlock = new block<T>[nrowblks * ncolblks];
    matrixBlock = new block[nrowblks * ncolblks];
    top = new int*[nrowblks];
    //top = (int **) malloc(nrowblks * sizeof(int *));
    nnzPerRow = (int *) malloc(nrowblks * sizeof(int));

    for(i = 0 ; i < nrowblks ; i++)
    {
        top[i] = new int[ncolblks];
        //top[i] = (int *) malloc(ncolblks * sizeof(int));
        nnzPerRow[i] = 0;
    }

    for(blkr = 0 ; blkr < nrowblks ; blkr++)
    {
        for(blkc = 0 ; blkc < ncolblks ; blkc++)
        {
            top[blkr][blkc] = 0;
            matrixBlock[blkr * ncolblks + blkc].nnz = 0;
        }
    }
    cout<<"finished memory allocation for block.."<<endl;


    //cout<<"calculatig nnz per block"<<endl;

    //calculatig nnz per block
    for(c = 0 ; c < numcols ; c++)
    {
        k1 = colptrs[c]+1;
        k2 = colptrs[c + 1] - 1+1;
        blkc = ceil((c + 1) / (double)wblk);
        //k1 = colptrs[c];
        //k2 = colptrs[c + 1] - 1;
        //blkc = ceil((c + 1) / (float)wblk);
        //cout<<"K1: "<<k1<<" K2: "<<k2<<" blkc: "<<blkc<<endl;

        for(k = k1 - 1 ; k < k2 ; k++)
        {
            r = irem[k]+1;
            //r = irem[k];
            blkr = ceil(r/(double)wblk);
            if((blkr - 1) >= nrowblks || (blkc - 1) >= ncolblks)
            {
                cout << "(" << blkr - 1 << ", " << blkc - 1 << ") doesn't exist" << endl;
            }
            else
            {
                matrixBlock[(blkr - 1) * ncolblks + (blkc - 1)].nnz++;  
            }    
        }
    }

    cout<<"finished counting nnz in each block"<<endl;

    for(blkc = 0 ; blkc < ncolblks; blkc++)
    {
        for(blkr = 0 ; blkr < nrowblks ; blkr++)
        {
            //cout<<"br: "<<blkr<<" bc: "<<blkc<<" roffset: "<<blkr*wblk<<" coffset: "<<blkc*wblk<<endl;
            matrixBlock[blkr * ncolblks + blkc].roffset = blkr * wblk + 1;
            matrixBlock[blkr * ncolblks + blkc].coffset = blkc * wblk + 1;
            //cout<<"here 1"<<endl;

            if(matrixBlock[blkr * ncolblks + blkc].nnz > 0)
            {
                nnzPerRow[blkr] += matrixBlock[blkr * ncolblks + blkc].nnz;
                //matrixBlock[blkr * ncolblks + blkc].rloc = new int[matrixBlock[blkr * ncolblks + blkc].nnz];
                matrixBlock[blkr * ncolblks + blkc].rloc = new unsigned short int[matrixBlock[blkr * ncolblks + blkc].nnz];
                //matrixBlock[blkr * ncolblks + blkc].cloc = new int[matrixBlock[blkr * ncolblks + blkc].nnz];
                matrixBlock[blkr * ncolblks + blkc].cloc = new unsigned short int[matrixBlock[blkr * ncolblks + blkc].nnz];
                //matrixBlock[blkr * ncolblks + blkc].val = new T[matrixBlock[blkr * ncolblks + blkc].nnz];
                matrixBlock[blkr * ncolblks + blkc].val = new double[matrixBlock[blkr * ncolblks + blkc].nnz];
            }
            else
            {
                matrixBlock[blkr * ncolblks + blkc].rloc = NULL;
                matrixBlock[blkr * ncolblks + blkc].cloc = NULL;
            }
        }
    }

    cout<<"allocating memory for each block"<<endl;

    //for(blkr=0;blkr<nrowblks;blkr++)
    //{
        //printf("nnzPerRow[%d] : %d\n", blkr, nnzPerRow[blkr]);
    //}
    //cout<<"end for"<<endl;

    printf("numrow = %d numcols = %d\n", numrows,numcols);

    for(c = 0 ; c < numcols ; c++)
    {
        k1 = colptrs[c]+1;
        k2 = colptrs[c + 1] - 1+1;
        blkc = ceil((c + 1) / (double)wblk);
        //k1 = colptrs[c];
        //k2 = colptrs[c + 1] - 1;
        //blkc = ceil((c + 1) / (float)wblk);

    //    printf("k1 = %d k2 = %d blkc = %d \n", k1,k2,blkc);

        for(k = k1 - 1 ; k < k2 ; k++)
        {
            r = irem[k]+1;
            //r = irem[k];

      //      printf("irem[%d] = %d\n",k , r);
            blkr = ceil(r / (double)wblk);

            matrixBlock[(blkr - 1) * ncolblks+blkc - 1].rloc[top[blkr-1][blkc-1]] = r - matrixBlock[(blkr - 1) * ncolblks + blkc - 1].roffset;
            matrixBlock[(blkr - 1) * ncolblks+blkc - 1].cloc[top[blkr-1][blkc-1]] = (c + 1) -  matrixBlock[(blkr - 1) * ncolblks + blkc - 1].coffset;
            matrixBlock[(blkr - 1) * ncolblks+blkc - 1].val[top[blkr-1][blkc-1]] = xrem[k];

            top[blkr-1][blkc-1]=top[blkr-1][blkc-1]+1;
        }
    }

    for(i = 0 ; i < nrowblks ; i++)
    {
        delete [] top[i];
    }
    delete [] top;
}

void splitProcs_LBCM(){

	int colklass = 0 ; 
	int rowklass = 0 ; 
	int diagklass = 0 ; 
	idiag = -1 ; 
	int kk = 0 ; 

	int mycol , myrow ; 

	int jj , rowkey ;

	int icol ; 

	for(icol = 0 ; icol < ndiag  ; icol++){
		jj = 0 ; 

		if(world_rank == kk ){

			idiag = icol + 1 ; 
			mycol = idiag ; 
			myrow = idiag ; 
			colklass = mycol ; 
			rowklass = myrow ; 
			diagklass = 1 ; 
		}

		kk++ ; 

		for(jj = 1 ; jj < nsegments ; jj++){

			if(world_rank == kk){

				mycol = icol + 1 ; 
				myrow = world_rank/nsegments  + jj ; 
				myrow = myrow % ndiag + 1 ; 
				colklass = mycol ;
				rowklass = myrow ; 

			}
			kk++ ; 
		}

	}


	MPI_Comm_split(MPI_COMM_WORLD , colklass , world_rank , &col_comm );
	MPI_Comm_rank(col_comm , &mysegment) ; 


	rowkey = world_rank % nsegments ; 
	MPI_Comm_split(MPI_COMM_WORLD , rowklass , rowkey , &row_comm);

	MPI_Comm_split(MPI_COMM_WORLD , diagklass , idiag , &diag_comm) ; 

}




void splitProcs(){


	ndiag = floor(sqrt(2.0*world_size)) ; 

	nsegments = (ndiag+1) / 2 ; 

	int err ; 

	if((2*nsegments) != (ndiag+1)){
		MPI_Abort(MPI_COMM_WORLD , err) ; 
	}

	npe = ndiag * (ndiag+1) ; 
	npe = npe / 2 ; 

	if(npe != world_size ){
		MPI_Abort(MPI_COMM_WORLD , err) ; 
	}

	splitProcs_LBCM(); 

}





int main(int argc , char *argv[]){

	int blocksize ; 
	int rhs ; 


	//stringstream s(argv[1]);
    //s >> rhs;
    //stringstream s1(argv[2]);
    //s1 >> blocksize;

    double *xrem;
    block *matrixBlock;

    wblk = blocksize ;

    //char *filename = argv[3] ; 

    //read_custom(filename , xrem) ; 
    //csc2blkcoord(matrixBlock , xrem);


/////////////////MPI Initialization ////////////////

    int provided ; 
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
//	int world_size ;

	MPI_Comm_size(MPI_COMM_WORLD,&world_size);

//	int world_rank ; 
	MPI_Comm_rank(MPI_COMM_WORLD , &world_rank ) ;





    ///////////////////////////create MPI communicators//////////

	splitProcs() ; 

	/////////////////////////////////////////////////////////////

	printf("myrank %d mysegment %d\n",world_rank , mysegment);


	if(world_rank == 0 ){

		printf("subvectors = %d\n",ndiag);
		printf("subvector segments = %d\n",nsegments);
		printf("total segment = %d\n", ndiag*nsegments);


	}	



	////////////allocate ilvecshare and ilvecdisp/////////////


	ilvecshare = (int*) malloc(nsegments * sizeof(int)) ; 
	ilvecdisp = (int*) malloc((nsegments+1) * sizeof(int)) ; 


	printf("main rank %d idiag %d \n", world_rank , idiag);



	////////open diag vector file displacement//////////


	char filename[100] ; 

	FILE *fp ; 


	if(idiag > 0 ){

		strcpy(filename,"node006/veshare00");
		char idiag_char = idiag+48; 
		strncat(filename , &idiag_char , 1) ; 
		fp = fopen(filename , "r") ; 

		printf("file %s opened\n", filename);

		for(int i = 0 ; i < nsegments ; i++){

			fscanf(fp,"%d",&ilvecshare[i]);
			printf("ilvecshare[%d] = %d\n", i , ilvecshare[i]);


		}


	}


	MPI_Bcast(ilvecshare , nsegments , MPI_INT , 0 , col_comm);

	ilvecdisp[0] = 0 ; 

	for(int i = 1 ; i <= nsegments ; i++){

		ilvecdisp[i] = ilvecshare[i-1] + ilvecdisp[i-1] ; 
	}

	printf("world_rank %d ilvecdisp %d ilvecshare %d\n", world_rank , ilvecshare[nsegments] , ilvecdisp[nsegments]);

	ncol = ilvecdisp[nsegments] + ilvecshare[nsegments] ; 

	lvecdim = ilvecshare[mysegment] ; 


	if(idiag > 0 ){

		printf("idiag %d ncol %d lvecdim %d\n",idiag , ncol , lvecdim) ; 

		nrow = ncol ; 
		totaldim = ncol ; 

		MPI_Allreduce(MPI_IN_PLACE , totaldim , 1 , MPI_INT , MPI_SUM , diag_comm) ; 


	}



	///////bcast along columns here ///////////////


	///// no direct bcast along columns 



	///////////////////////////////////////////////




	///////////// bcast along rows here /////////////////

	MPI_Bcast(nrow , 1 , MPI_INT , 0 , row_comm) ; 

	////////////////////////////////////////////////////





/////////////////// SPMM ////////////////////////////////

	#pragma omp parallel for default(shared) private(tstart, rbase, cbase, j, k, l, r, c, xcoef, index)
    for(i = 0; i < nrowblks ; i++)
    {
        //tstart = omp_get_wtime();
        rbase = H[i * ncolblks + 0].roffset;
        
        for(j = 0 ; j < ncolblks ; j++)
        {
            index = i * ncolblks + j;
            //cbase = H[i * ncolblks + j].coffset;
            cbase = H[index].coffset;
            //if(H[i * ncolblks + j].nnz > 0)
            if(H[index].nnz > 0)
            {
                //for(k = 0 ; k < H[i * ncolblks + j].nnz ; k++)
                for(k = 0 ; k < H[index].nnz ; k++)
                {
                    //r = rbase + H[i * ncolblks + j].rloc[k] - 1;
                    r = ( rbase + H[index].rloc[k] - 1 ) << 6;
                    //c = cbase + H[i * ncolblks + j].cloc[k] - 1;
                    c = ( cbase + H[index].cloc[k] - 1 ) << 6;
                    //xcoef = H[i * ncolblks + j].val[k];
                    xcoef = H[index].val[k];
                    #pragma omp simd 
                    //for(l = 0; l < blocksize ; l++)
                    for(l = 0; l < 8 ; l++)
                    {
                        //Y[r * blocksize + l] = Y[r * blocksize + l] + xcoef * X[c * blocksize + l];
                        Hammp[r + l] = Hammp[r + l] + xcoef * amp[c + l];
                    }
                }
            }
        }
        //taskTiming[omp_get_thread_num()][1] += (omp_get_wtime() - tstart);
    } //end for



//////////////////////////// Reduce the results along row communicators ////////////////////////

	if(idiag > 0){


		MPI_Reduce(MPI_IN_PLACE , Hammp , 8 * nrow , MPI_DOUBLE , MPI_SUM , 0 , row_comm) ; 



	}

	else{


		MPI_Reduce(Hammp , 0 , 8*nrow , MPI_DOUBLE , MPI_SUM , 0 , row_comm) ; 

	}




    /////////////////////////////////////////////////////////////////////////////////////





    ///////////////////////// SPMMT //////////////////////////////////


    #pragma omp parallel for default(shared) private(tstart, rbase, cbase, j, k, l, r, c, xcoef, index)
    for(j = 0; j < ncolblks ; j++)
    {
        //tstart = omp_get_wtime();
        cbase = H[j * nrowblks + 0].roffset;
        
        for(i = 0 ; i < nrowblks ; i++)
        {
            index = j * nrowblks + i;
            //cbase = H[i * ncolblks + j].coffset;
            rbase = H[index].coffset;
            //if(H[i * ncolblks + j].nnz > 0)
            if(H[index].nnz > 0)
            {
                //for(k = 0 ; k < H[i * ncolblks + j].nnz ; k++)
                for(k = 0 ; k < H[index].nnz ; k++)
                {
                    //r = rbase + H[i * ncolblks + j].rloc[k] - 1;
                    r = ( rbase + H[index].rloc[k] - 1 ) << 6;
                    //c = cbase + H[i * ncolblks + j].cloc[k] - 1;
                    c = ( cbase + H[index].cloc[k] - 1 ) << 6;
                    //xcoef = H[i * ncolblks + j].val[k];
                    xcoef = H[index].val[k];
                    #pragma omp simd 
                    //for(l = 0; l < blocksize ; l++)
                    for(l = 0; l < 8 ; l++)
                    {
                        //Y[r * blocksize + l] = Y[r * blocksize + l] + xcoef * X[c * blocksize + l];
                       	HammpT[c + l] = HammpT[c + l] + xcoef * ampT[r + l];
                    }
                }
            }
        }
        //taskTiming[omp_get_thread_num()][1] += (omp_get_wtime() - tstart);
    } //end for



//////////////////////////////////////////



    /////////////////////////// reduce scatter //////////////////////////////


	MPI_Reduce_Scatter(HammpT , distHamp , nblk * ilvecshare , MPI_DOUBLE , MPI_SUM , col_comm) ; 



    //////////////////////////////////////////////////////////////////////////



	MPI_Finalize();













}