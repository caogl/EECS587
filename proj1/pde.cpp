#include <mpi.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <math.h>
using namespace std;


void printMatrix(vector<vector<double> >& matrix);

// return the initialized matrix with specific values
vector<vector<double> > initMatrix(int n, int p, int rank);
// calculate the mean of the neighbour
double meanNeighbour(vector<vector<double> >& matrix, int i, int j);

int main(int argc, char **argv)
{
	if(argc != 2)
	{
		cout << "Error, please specify the size of the matrix" << endl;
		exit(1);
	}
	
	int rank;
	int p;// number of processors
	int n=atoi(argv[1]);

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	// init the matrix with given values
	vector<vector<double> > matrix=initMatrix(n, p, rank);

	MPI_Barrier(MPI_COMM_WORLD);
	double startTime=MPI_Wtime();
	

	
	// begin  the 500 iterations
	for(int k=0; k<500; k++)
	{

		// send the updated value rightwards and receive leftwards
                if(rank!=0)
                {
                	MPI_Recv(&matrix[0].front(), n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
                MPI_Send(&matrix[matrix.size()-2].front(), n, MPI_DOUBLE, (rank+1)%p, 0, MPI_COMM_WORLD);
		
                if(rank==0)
                {
                        MPI_Recv(&matrix[0].front(), n, MPI_DOUBLE, p-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

		// send the updated value leftwards and receive rightwards
		if(rank!=p-1)
		{
			MPI_Recv(&matrix[matrix.size()-1].front(), n, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if(rank!=0)
			MPI_Send(&matrix[1].front(), n, MPI_DOUBLE, (rank-1), 1, MPI_COMM_WORLD);
		else
			MPI_Send(&matrix[1].front(), n, MPI_DOUBLE, p-1, 1, MPI_COMM_WORLD);
		if(rank==p-1)
		{
			MPI_Recv(&matrix[0].front(), n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
			
		// Doing the mean-neighbour computation
		for(int i=1; i<matrix.size()-1; i++)
		{
			for(int j=1; j<matrix[i].size()-1; j++)
			{
				matrix[i][j]=meanNeighbour(matrix, i, j);
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double endTime, totalTime;
	if(rank==0)
	{
		double endTime=MPI_Wtime();
		double totalTime = endTime - startTime;
		cout<<"The total time is: "<<totalTime<<endl;
	}
	
	/*************verify the sum using reduce***************/
	// calculate local sum for each processor
	int bandwidth=n/p;
	double localSum;
	double globalSum;
	for(int i=1; i<matrix.size()-1; i++)
	{
		localSum+=matrix[i][bandwidth*rank+i-1];
	}
	MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank==0)
	{
		cout<<"The diagonal sum is: "<<globalSum<<endl;
		//cout<<"The total run time is: "<<totalTime<<endl;
	}


	MPI_Finalize();
	return 0;
	
}

vector<vector<double> > initMatrix(int n, int p, int rank)
{
        int bandwidth=n/p;
        int bandwidthActual=bandwidth;
        if(rank==p-1)
                bandwidthActual=n-(p-1)*bandwidth;

        vector<double> tmp(n, 0.5);
        vector<vector<double> > result(bandwidthActual, tmp);
        double temp;

	int j=0;
        for(int i=rank*bandwidth; i<rank*bandwidth+bandwidthActual; i++)
        {
                result[j][0]=0;
                temp=pow((double)i/n,2);
                temp=5*sin(temp*M_PI);
                result[j][n-1]=temp;
		j++;
        }

	result.insert(result.begin(), tmp);
	result.push_back(tmp);
	
        return result;
}

double meanNeighbour(vector<vector<double> >& matrix, int i, int j)
{
	double sum=matrix[i][j];
	sum+=matrix[i-1][j];
        sum+=matrix[i+1][j];
        sum+=matrix[i-1][j-1];
        sum+=matrix[i][j-1];
        sum+=matrix[i+1][j-1];
        sum+=matrix[i-1][j+1];
        sum+=matrix[i][j+1];
        sum+=matrix[i+1][j+1];
	return sum/9;	
}

void printMatrix(vector<vector<double> >& matrix)
{
	cout<<"begin to print"<<endl;
	cout<<"Size of the matrix is: "<<matrix.size()<<" "<<matrix[0].size()<<endl;
	for(int i=0; i<matrix.size(); i++)
	{
		for(int j=0; j<matrix[i].size(); j++)
			cout<<matrix[i][j]<<" ";
		cout<<endl;
	}

	cout<<endl;
}









