#define ROWS  4
#define COLUMNS 4
#define MASTER 0

typedef struct {
	int id;
	float height;
	float vol;

} cuboid;

void getResultsAndWriteToFile(int size, int id, int lineAndColSize,
		cuboid *myCuboid, MPI_Datatype *cuboidTransfer, MPI_Comm *comm);
void readFromFile(cuboid **cuboidArray, int size, char **fileName);
void printCuboidArray(int size, cuboid *cuboidArray);
void lineTransaction(cuboid *myCuboid, int *coord, int dest, int source,
		MPI_Datatype *cuboidTransfer, MPI_Comm *comm, MPI_Status *status,
		int id, void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
		void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid));
void colTransaction(cuboid *myCuboid, int *coord, int dest, int source,
		MPI_Datatype *cuboidTransfer, MPI_Comm *comm, MPI_Status *status,
		int id, void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
		void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid));
void oddEven(int n, cuboid *myCuboid, int *coord, int dest, int source,
		MPI_Datatype *cuboidTransfer, MPI_Comm *comm, MPI_Status *status,
		void (*transFunc)(cuboid *myCuboid, int *coord, int dest, int source,
				MPI_Datatype *cuboidTransfer, MPI_Comm *comm,
				MPI_Status *status, int id,
				void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
				void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid)),
		int id);
void getMinCuboid(cuboid *myCuboid, cuboid *tempCuboid);
void getMaxCuboid(cuboid *myCuboid, cuboid *tempCuboid);
void replaceMyCuboid(cuboid *myCuboid, cuboid *tempCuboid);
void init(int *argc, char **argv, int *rank, int *size, int *n);
void createDataType(cuboid myCuboid, MPI_Datatype *cuboidTransfer);
void createCartesian(int lineAndRowSize, MPI_Comm *comm);
void writeResultToFile(int n, cuboid *result);
void allocateCuboidArray(int size, cuboid **result);
void mpiShearSortPreparations(int rank, int *coord, int *id, int *source,
		int *dest, char **argv, cuboid *myCuboid, MPI_Datatype *cuboidTransfer,
		MPI_Comm *comm, cuboid *cuboidArray);

void shearSort(int n, int id, cuboid *myCuboid, int *coord, int *dest,
		int *source, MPI_Datatype *cuboidTransfer, MPI_Comm *comm,
		MPI_Status *status,
		void (**transFunc)(cuboid *myCuboid, int *coord, int dest, int source,
				MPI_Datatype *cuboidTransfer, MPI_Comm *comm,
				MPI_Status *status, int id,
				void (*compareAndSwap_1)(cuboid *myCuboid, cuboid *tempCuboid),
				void (*compareAndSwap_2)(cuboid *myCuboid, cuboid *tempCuboid)));
int properExit();
