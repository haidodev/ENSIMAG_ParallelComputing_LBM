/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
//
// GOAL: Implement a 1D communication scheme along
//       X axis with blocking communications.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

enum {
	TAG_LEFT,
	TAG_RIGHT,
};

/****************************************************/
void lbm_comm_init_ex1(lbm_comm_t * comm, int total_width, int total_height)
{
	//get infos
	int rank;
	int comm_size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

	// check
	if (total_width % comm_size != 0)
		fatal("Invalid total width, should be a multiple of the number of processes !");

	comm->nb_x = comm_size;
	comm->nb_y = 1;

	comm->rank_x = rank;
	comm->rank_y = 0;

	comm->width = (total_width / comm_size) + 2;
	comm->height = total_height + 2;
	printf("rank %d : width = %d, height = %d\n", rank, comm->width, comm->height);
	// TODO : calculate the absolute position in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = rank * (total_width / comm_size);
	comm->y = 0;

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex1(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with blocking MPI functions (MPI_Send & MPI_Recv)
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[DIRECTIONS] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	if (comm->rank_x > 0)
	{
		// MPI_Sendrecv(
		// 	lbm_mesh_get_cell(mesh, 1, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x - 1, 0,
		// 	lbm_mesh_get_cell(mesh, 0, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x - 1, 0,
		// 	MPI_COMM_WORLD, MPI_STATUS_IGNORE
		// );

		// This implementation is may block if the height is too large (detail in the report), you can use the following implementation to avoid this issue (but less efficient) :
		MPI_Send(lbm_mesh_get_cell(mesh, 1, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x - 1, TAG_LEFT, MPI_COMM_WORLD);
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x - 1, TAG_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// This implementation is suboptimal since we can send all the data at once, reduce the latency
		// for (int i = 0; i < comm->height; i++)
		// {
		// 	MPI_Send(lbm_mesh_get_cell(mesh, 1, i), DIRECTIONS, MPI_DOUBLE, comm->rank_x - 1, TAG_LEFT, MPI_COMM_WORLD);
		// 	MPI_Recv(lbm_mesh_get_cell(mesh, 0, i), DIRECTIONS, MPI_DOUBLE, comm->rank_x - 1, TAG_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// }
	}
	if (comm->rank_x < comm->nb_x - 1)
	{
		// MPI_Sendrecv(
		// 	lbm_mesh_get_cell(mesh, comm->width - 2, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x + 1, 0,
		// 	lbm_mesh_get_cell(mesh, comm->width - 1, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x + 1, 0,
		// 	MPI_COMM_WORLD, MPI_STATUS_IGNORE
		// );
		MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x + 1, TAG_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, 0), DIRECTIONS * comm->height, MPI_DOUBLE, comm->rank_x + 1, TAG_RIGHT, MPI_COMM_WORLD);
		// for (int i = 0; i < comm->height; i++)
		// {
		// 	MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, i), DIRECTIONS, MPI_DOUBLE, comm->rank_x + 1, TAG_RIGHT, MPI_COMM_WORLD);
		// 	MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, i), DIRECTIONS, MPI_DOUBLE, comm->rank_x + 1, TAG_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// }
	}

}
