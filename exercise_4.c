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
//       8 neighbors using manual copy for non
//       contiguous side and blocking communications
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
//     - Manual copy for non continguous cells
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"
#include <math.h>
#include <string.h>

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
void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
	// //get infos
	// int rank;
	// int comm_size;
	// MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	// MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

	// // check
	// if (total_width % comm_size != 0)
	// 	fatal("Invalid total width, should be a multiple of the number of processes !");

	// comm->nb_x = comm_size;
	// comm->nb_y = 1;

	// comm->rank_x = rank;
	// comm->rank_y = 0;

	// comm->width = (total_width / comm_size) + 2;
	// comm->height = total_height + 2;
	// printf("rank %d : width = %d, height = %d\n", rank, comm->width, comm->height);
	// // TODO : calculate the absolute position in the global mesh.
	// //        without accounting the ghost cells
	// //        (used to setup the obstable & initial conditions).
	// comm->x = rank * (total_width / comm_size);
	// comm->y = 0;

	// //if debug print comm
	// //lbm_comm_print(comm);

	//get infos
	int rank;
	int comm_size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

	//
	// TODO: calculate the splitting parameters for the current task.
	//
	
	// assume square domain and square splitting for simplicity
	int nb_tasks_per_dim = (int) sqrt(comm_size);
	if (nb_tasks_per_dim * nb_tasks_per_dim != comm_size)
		fatal("Invalid number of processes, should be a perfect square !");
	if (total_width % nb_tasks_per_dim != 0 || total_height % nb_tasks_per_dim != 0)
		fatal("Invalid total width/height, should be a multiple of the number of tasks per dimension !");

	// TODO: calculate the number of tasks along X axis and Y axis.
	comm->nb_x = nb_tasks_per_dim;
	comm->nb_y = nb_tasks_per_dim;

	// TODO: calculate the current task position in the splitting
	comm->rank_x = rank % nb_tasks_per_dim;
	comm->rank_y = rank / nb_tasks_per_dim;

	// TODO : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = (total_width / nb_tasks_per_dim) + 2;
	comm->height = (total_height / nb_tasks_per_dim) + 2;

	// TODO : calculate the absolute position  (in cell number) in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = comm->rank_x * (total_width / nb_tasks_per_dim);
	comm->y = comm->rank_y * (total_height / nb_tasks_per_dim);

	//OPTIONAL : if you want to avoid allocating temporary copy buffer
	//           for every step :
	//comm->buffer_recv_down, comm->buffer_recv_up, comm->buffer_send_down, comm->buffer_send_up
	comm -> buffer_send_up = (double *) malloc(sizeof(double) * DIRECTIONS * comm->width);
	comm -> buffer_recv_up = (double *) malloc(sizeof(double) * DIRECTIONS * comm->width);
	comm -> buffer_send_down = (double *) malloc(sizeof(double) * DIRECTIONS * comm->width);
	comm -> buffer_recv_down = (double *) malloc(sizeof(double) * DIRECTIONS * comm->width);

	//if debug print comm
	// lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
	//free allocated ressources
	free(comm -> buffer_send_up);
	free(comm -> buffer_recv_up);
	free(comm -> buffer_send_down);
	free(comm -> buffer_recv_down);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - manual copy in temp buffer for non contiguous side 
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

	// Communication with left and right neighbors (contiguous)
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

	if (comm->rank_y > 0)
	{
		int start = comm -> rank_x > 0 ? 1 : 0, end = comm -> rank_x < comm -> nb_x - 1 ? comm -> width - 1 : comm -> width;
		for (int i = start; i < end; i++)
			memcpy(&(comm -> buffer_send_up[i * DIRECTIONS]), lbm_mesh_get_cell(mesh, i, 1), sizeof(double) * DIRECTIONS);

		int up = (comm->rank_y - 1) * comm->nb_x + comm->rank_x;
		MPI_Send(comm->buffer_send_up, DIRECTIONS * comm->width, MPI_DOUBLE, up, TAG_UP, MPI_COMM_WORLD);
		MPI_Recv(comm->buffer_recv_up, DIRECTIONS * comm->width, MPI_DOUBLE, up, TAG_DOWN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int i = start; i < end; i++)
			memcpy(lbm_mesh_get_cell(mesh, i, 0), &(comm -> buffer_recv_up[i * DIRECTIONS]), sizeof(double) * DIRECTIONS);
	}

	if (comm->rank_y < comm->nb_y - 1)
	{
		int start = comm -> rank_x > 0 ? 1 : 0, end = comm -> rank_x < comm -> nb_x - 1 ? comm -> width - 1 : comm -> width;
		for (int i = start; i < end; i++)
			memcpy(&(comm -> buffer_send_down[i * DIRECTIONS]), lbm_mesh_get_cell(mesh, i, comm->height - 2), sizeof(double) * DIRECTIONS);

		int down = (comm->rank_y + 1) * comm->nb_x + comm->rank_x;
		MPI_Recv(comm->buffer_recv_down, DIRECTIONS * comm->width, MPI_DOUBLE, down, TAG_UP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(comm->buffer_send_down, DIRECTIONS * comm->width, MPI_DOUBLE, down, TAG_DOWN, MPI_COMM_WORLD);
		for (int i = start; i < end; i++)
			memcpy(lbm_mesh_get_cell(mesh, i, comm->height - 1), &(comm -> buffer_recv_down[i * DIRECTIONS]), sizeof(double) * DIRECTIONS);
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
