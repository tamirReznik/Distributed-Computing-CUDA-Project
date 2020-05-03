//Made by: Tamir Reznik

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "HW5.h"

int main(int argc, char *argv[]) {
	//Array of 2 odd even sorting functions - [0] - line , [1] - col
	void (*transitionFunction[2])(cuboid *myCuboid, int *coord, int dest,
			int source, MPI_Datatype *cuboidTransfer, MPI_Comm *comm,
			MPI_Status *status, int id,
			void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
			void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid)) = {lineTransaction, colTransaction
	};
	int rank, size, id, lineAndColSize;
	int coord[2], source[2], dest[2];
	MPI_Comm comm;
	MPI_Status status;
	cuboid myCuboid;
	MPI_Datatype cuboidTransfer;

	init(&argc, argv, &rank, &size, &lineAndColSize);
	cuboid cuboidArray[size];

// Create MPI user data type for cuboid
	createDataType(myCuboid, &cuboidTransfer);

// A two-dimensional cylinder of 'size' processes in a square_root(size) grid
	createCartesian(&lineAndColSize, &comm);

//	Only master read from file
	if (rank == MASTER)
		readFromFile(cuboidArray, size, argv);

	mpiShearSortPreparations(&rank, coord, &id, source, dest, argv, &myCuboid,
			&cuboidTransfer, &comm, cuboidArray);

	shearSort(&lineAndColSize, &id, &myCuboid, coord, dest, source,
			&cuboidTransfer, &comm, &status, transitionFunction);

	getResultsAndWriteToFile(&size, &id, &lineAndColSize, &myCuboid,
			&cuboidTransfer, &comm);

	return properExit();
}

void getResultsAndWriteToFile(int *size, int *id, int *lineAndColSize,
		cuboid *myCuboid, MPI_Datatype *cuboidTransfer, MPI_Comm *comm) {
	cuboid *result;

	allocateResultArray(*size, &result);

	MPI_Gather(myCuboid, 1, *cuboidTransfer, result, 1, *cuboidTransfer, MASTER,
			*comm);

	if (*id == MASTER)
		writeResultToFile(*lineAndColSize, result);

	free(result);
}

int properExit() {
	fflush(stdout);
	MPI_Finalize();
	return 0;
}

void shearSort(int *n, int *id, cuboid *myCuboid, int *coord, int *dest,
		int *source, MPI_Datatype *cuboidTransfer, MPI_Comm *comm,
		MPI_Status *status,
		void (**transFunc)(cuboid *myCuboid, int *coord, int dest, int source,
				MPI_Datatype *cuboidTransfer, MPI_Comm *comm,
				MPI_Status *status, int id,
				void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
				void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid))) {
	int i;
	int numOfIteration = (int) (2 * log(*n) / log(2) + 1);
	for (i = 0; i < numOfIteration; ++i)
		oddEven(*n, myCuboid, coord, dest[i % 2], source[i % 2], cuboidTransfer,
				comm, status, transFunc[i % 2], *id);
}

void mpiShearSortPreparations(int *rank, int *coord, int *id, int *source,
		int *dest, char **argv, cuboid *myCuboid, MPI_Datatype *cuboidTransfer,
		MPI_Comm *comm, cuboid *cuboidArray) {

	MPI_Scatter(cuboidArray, 1, *cuboidTransfer, myCuboid, 1, *cuboidTransfer,
	MASTER,
	MPI_COMM_WORLD);

	// Each process displays its rank and cartesian coordinates
	MPI_Cart_coords(*comm, *rank, 2, coord);

	MPI_Cart_rank(*comm, coord, id);

	MPI_Cart_shift(*comm, 1, 1, &source[1], &dest[1]);
	if (coord[1] % 2 == 0)
		MPI_Cart_shift(*comm, 0, 1, &source[0], &dest[0]);
	else
		MPI_Cart_shift(*comm, 0, 1, &dest[0], &source[0]);

}

void allocateResultArray(int size, cuboid **result) {
	*result = (cuboid*) malloc(size * sizeof(cuboid));
	if (*result == NULL) {
		printf("Memory not allocated.\n");
		exit(0);
	}
}

void writeResultToFile(int n, cuboid *result) {
	int i, j;

	FILE *output = fopen("result.dat", "w+t");
	if (output == NULL) {
		printf("Cannot open result.dat file...\n");
		exit(1);
	}
	//write to file...
	for (i = 0; i < n; ++i)
		for (j = 0; j < n; ++j)
			if (i % 2 == 0)
				fprintf(output, "%d  ", result[j * n + i].id);
			else
				fprintf(output, "%d  ", result[(n - j - 1) * n + i].id);

	if (fclose(output) != 0) {
		printf("Cannot close result.dat file...\n");
		exit(1);
	}

}

//Create Cartesian mpi group via mpi library
void createCartesian(int *n, MPI_Comm *comm) {
	int dim[2], period[2] = { 0 }, reorder;
	dim[0] = *n;
	dim[1] = *n;
	period[0] = 0;
	period[1] = 0;
	reorder = 0;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, comm);
}

void init(int *argc, char **argv, int *rank, int *size, int *n) {
	MPI_Init(argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, rank);
	MPI_Comm_size(MPI_COMM_WORLD, size);
	*n = sqrt(*size);
//	*numOfIteration = (int) (2 * log(*n) / log(2) + 1);

	if (*size != (*n) * (*n) || *argc != 2) {

		printf(
				"Please enter a number whose root is an integer\nMake sure only to call the program with file name as first argument");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

void createDataType(cuboid myCuboid, MPI_Datatype *cuboidTransfer) {

	MPI_Aint disp[3];
	MPI_Datatype type[3] = { MPI_INT, MPI_FLOAT, MPI_FLOAT };
	int blocklen[3] = { 1, 1, 1 };
	// Create MPI user data type for cuboid
	disp[0] = (char*) &myCuboid.id - (char*) &myCuboid;
	disp[1] = (char*) &myCuboid.height - (char*) &myCuboid;
	disp[2] = (char*) &myCuboid.vol - (char*) &myCuboid;
	MPI_Type_create_struct(3, blocklen, disp, type, cuboidTransfer);
	MPI_Type_commit(cuboidTransfer);

}

void readFromFile(cuboid *cuboidArray, int size, char **fileName) {
	int i;
	float width, length;

	FILE *cuboidFile = fopen(fileName[1], "r");
	if (cuboidFile == NULL) {
		fprintf(stderr, "line %d -->Cannot open file Try again later.\n",
		__LINE__);
		perror("fopen ");
		exit(1);
	}
	for (i = 0; i < size; ++i) {
		if (fscanf(cuboidFile, "%d%f%f%f", &cuboidArray[i].id, &length, &width,
				&cuboidArray[i].height) != 4) {
			fprintf(stderr,
					"line %d -->Cannot read specific cuboid from file.\n",
					__LINE__);
			exit(1);
		};
		cuboidArray[i].vol = length * cuboidArray[i].height * width;
	}
}

void printCuboidArray(int size, cuboid *cuboidArray) {
	int i;
	for (i = 0; i < size; ++i)
		printf("id: %d, height: %f vol: %f\n", cuboidArray[i].id,
				cuboidArray[i].height, cuboidArray[i].vol);
}

//set myCuboid to be the minimum cuboid between the two
void getMinCuboid(cuboid *myCuboid, cuboid *tempCuboid) {

	if (myCuboid->vol > tempCuboid->vol)
		replaceMyCuboid(myCuboid, tempCuboid);

	else if (myCuboid->vol == tempCuboid->vol)

		if (myCuboid->height > tempCuboid->height)
			replaceMyCuboid(myCuboid, tempCuboid);

}

//set myCuboid to be the maximum cuboid between the two
void getMaxCuboid(cuboid *myCuboid, cuboid *tempCuboid) {
	if (myCuboid->vol < tempCuboid->vol)
		replaceMyCuboid(myCuboid, tempCuboid);

	else if (myCuboid->vol == tempCuboid->vol)

		if (myCuboid->height < tempCuboid->height)
			replaceMyCuboid(myCuboid, tempCuboid);
}

void replaceMyCuboid(cuboid *myCuboid, cuboid *tempCuboid) {
	myCuboid->id = tempCuboid->id;
	myCuboid->height = tempCuboid->height;
	myCuboid->vol = tempCuboid->vol;
}

//odd-Even for column in the matrix
void colTransaction(cuboid *myCuboid, int *coord, int dest, int source,
		MPI_Datatype *cuboidTransfer, MPI_Comm *comm, MPI_Status *status,
		int id, void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
		void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid)) {

	cuboid tempCuboid;
	if (coord[1] % 2 == 0) {
		if (dest > -1) {
			MPI_Send(myCuboid, 1, *cuboidTransfer, dest, 0, *comm);

			MPI_Recv(&tempCuboid, 1, *cuboidTransfer, dest, 0, *comm, status);

			(*compareAndSwap_1)(myCuboid, &tempCuboid);
		}
	} else if (source > -1) {
		MPI_Recv(&tempCuboid, 1, *cuboidTransfer, source, 0, *comm, status);

		MPI_Send(myCuboid, 1, *cuboidTransfer, source, 0, *comm);

		(*compareAndSwap_2)(myCuboid, &tempCuboid);

	}

}

//odd-Even for line in the matrix
void lineTransaction(cuboid *myCuboid, int *coord, int dest, int source,
		MPI_Datatype *cuboidTransfer, MPI_Comm *comm, MPI_Status *status,
		int id, void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
		void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid)) {
	cuboid tempCuboid;
	if (coord[0] % 2 == coord[1] % 2) {
		if (dest > -1) {
			MPI_Send(myCuboid, 1, *cuboidTransfer, dest, 0, *comm);

			MPI_Recv(&tempCuboid, 1, *cuboidTransfer, dest, 0, *comm, status);

			(*compareAndSwap_1)(myCuboid, &tempCuboid);

		}
	} else if (source > -1) {
		MPI_Recv(&tempCuboid, 1, *cuboidTransfer, source, 0, *comm, status);

		MPI_Send(myCuboid, 1, *cuboidTransfer, source, 0, *comm);

		(*compareAndSwap_2)(myCuboid, &tempCuboid);

	}

}

//generic odd-Even for cold and rows
void oddEven(int n, cuboid *myCuboid, int *coord, int dest, int source,
		MPI_Datatype *cuboidTransfer, MPI_Comm *comm, MPI_Status *status,
		void (*transFunc)(cuboid *myCuboid, int *coord, int dest, int source,
				MPI_Datatype *cuboidTransfer, MPI_Comm *comm,
				MPI_Status *status, int id,
				void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
				void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid)),
		int id) {
	int j;
	for (j = 0; j < n; ++j)
		if (j % 2 == 0)

			(*transFunc)(myCuboid, coord, dest, source, cuboidTransfer, comm,
					status, id, getMinCuboid, getMaxCuboid);

		else

			(*transFunc)(myCuboid, coord, source, dest, cuboidTransfer, comm,
					status, id, getMaxCuboid, getMinCuboid);

}

