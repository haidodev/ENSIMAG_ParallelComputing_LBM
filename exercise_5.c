/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication scheme with
//      8 neighbors using MPI types for non contiguous
//      side.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
// NEW:
//     - >>> MPI type for non contiguous cells <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

enum {
	TAG_UP,
	TAG_DOWN,
	TAG_LEFT,
	TAG_RIGHT,
	TAG_UP_LEFT,
	TAG_UP_RIGHT,
	TAG_DOWN_LEFT,
	TAG_DOWN_RIGHT
};
/****************************************************/
void lbm_comm_init_ex5(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5 execpt for type creation
	lbm_comm_init_ex4(comm, total_width, total_height);

	//TODO: create MPI type for non contiguous side in comm->type.
	MPI_Type_vector(comm -> width - 2, DIRECTIONS, comm -> height * DIRECTIONS, MPI_DOUBLE, &comm -> type);
	MPI_Type_commit(&comm -> type);
}

/****************************************************/
void lbm_comm_release_ex5(lbm_comm_t * comm)
{
	//we use the same implementation than ex5 except for type destroy
	lbm_comm_release_ex4(comm);

	//TODO: release MPI type created in init.
	MPI_Type_free(&comm -> type);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex5(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - use MPI type for non contiguous side 
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	//
	// TIP: create a function to get the target rank from x,y task coordinate.
	// TIP: You can use MPI_PROC_NULL on borders.
	// TIP: send the corner values 2 times, with the up/down/left/write communication
	//      and with the diagonal communication in a second time, this avoid
	//      special cases for border tasks.

	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications
	if (comm->rank_x > 0)
	{
		int start = comm -> rank_y > 0 ? 1 : 0, end = comm -> rank_y < comm -> nb_y - 1 ? comm -> height - 1 : comm -> height;
		int elements_to_send = (end - start) * DIRECTIONS;
		MPI_Send(lbm_mesh_get_cell(mesh, 1, start), elements_to_send, MPI_DOUBLE, comm -> rank_y * comm->nb_x + comm->rank_x - 1, TAG_LEFT, MPI_COMM_WORLD);
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, start), elements_to_send, MPI_DOUBLE, comm -> rank_y * comm->nb_x + comm->rank_x - 1, TAG_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if (comm->rank_x < comm->nb_x - 1)
	{
		int start = comm -> rank_y > 0 ? 1 : 0, end = comm -> rank_y < comm -> nb_y - 1 ? comm -> height - 1 : comm -> height;
		int elements_to_send = (end - start) * DIRECTIONS;
		MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, start), elements_to_send, MPI_DOUBLE, comm -> rank_y * comm->nb_x + comm->rank_x + 1, TAG_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, start), elements_to_send, MPI_DOUBLE, comm -> rank_y * comm->nb_x + comm->rank_x + 1, TAG_RIGHT, MPI_COMM_WORLD);
	}


	// Communication with top and bottom neighbors (non contiguous, need manual copy)
	if (comm->rank_y > 0)
	{
		MPI_Send(lbm_mesh_get_cell(mesh, 1, 1), 1, comm -> type, (comm->rank_y - 1) * comm->nb_x + comm->rank_x, TAG_UP, MPI_COMM_WORLD);
		MPI_Recv(lbm_mesh_get_cell(mesh, 1, 0), 1, comm -> type, (comm->rank_y - 1) * comm->nb_x + comm->rank_x, TAG_DOWN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (comm -> rank_x == 0)
		{
			MPI_Send(lbm_mesh_get_cell(mesh, 0, 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x, TAG_UP, MPI_COMM_WORLD);
			MPI_Recv(lbm_mesh_get_cell(mesh, 0, 0), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x, TAG_DOWN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		if (comm -> rank_x == comm -> nb_x - 1)
		{
			MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 1, 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x, TAG_UP, MPI_COMM_WORLD);
			MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, 0), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x, TAG_DOWN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

	}
	if (comm->rank_y < comm->nb_y - 1)
	{
		MPI_Recv(lbm_mesh_get_cell(mesh, 1, comm->height - 1), 1, comm -> type, (comm->rank_y + 1) * comm->nb_x + comm->rank_x, TAG_UP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(lbm_mesh_get_cell(mesh, 1, comm->height - 2), 1, comm -> type, (comm->rank_y + 1) * comm->nb_x + comm->rank_x, TAG_DOWN, MPI_COMM_WORLD);
		if (comm -> rank_x == 0)
		{
			MPI_Recv(lbm_mesh_get_cell(mesh, 0, comm->height - 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x, TAG_UP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(lbm_mesh_get_cell(mesh, 0, comm->height - 2), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x, TAG_DOWN, MPI_COMM_WORLD);
		}

		if (comm -> rank_x == comm -> nb_x - 1)
		{
			MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, comm->height - 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x, TAG_UP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 1, comm->height - 2), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x, TAG_DOWN, MPI_COMM_WORLD);
		}

	}

	// Communication with diagonal neighbors (contiguous)
	if (comm->rank_x > 0 && comm->rank_y > 0)
	{
		MPI_Send(lbm_mesh_get_cell(mesh, 1, 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x - 1, TAG_UP_LEFT, MPI_COMM_WORLD);
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, 0), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x - 1, TAG_DOWN_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if (comm->rank_x < comm->nb_x - 1 && comm->rank_y > 0)
	{
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x + 1, TAG_UP_RIGHT, MPI_COMM_WORLD);
		MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, 0), DIRECTIONS, MPI_DOUBLE, (comm->rank_y - 1) * comm->nb_x + comm->rank_x + 1, TAG_DOWN_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if (comm->rank_x > 0 && comm->rank_y < comm->nb_y - 1)
	{
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, comm->height - 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x - 1, TAG_UP_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(lbm_mesh_get_cell(mesh, 1, comm->height - 2), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x - 1, TAG_DOWN_LEFT, MPI_COMM_WORLD);
	}
	if (comm->rank_x < comm->nb_x - 1 && comm->rank_y < comm->nb_y - 1)
	{
		MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, comm->height - 1), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x + 1, TAG_UP_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, comm->height - 2), DIRECTIONS, MPI_DOUBLE, (comm->rank_y + 1) * comm->nb_x + comm->rank_x + 1, TAG_DOWN_RIGHT, MPI_COMM_WORLD);
	}

}
