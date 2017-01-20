#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#define DEFAULT_START_X 0.0
#define DEFAULT_END_X 2.0
#define DEFAULT_START_Y 0.0
#define DEFAULT_END_Y 2.0
#define DEFAULT_VALUE_P 2.5
#define DEFAULT_COEFF 0.66666
//1.5

// Объявление переменных
int global_length;
int global_width;
double x_coord;
double y_coord;
int x_bias;
int y_bias;
double xl_step;
double xr_step;
double yt_step;
double yb_step;
double exponential_coeff;
double* x_coords;
double* y_coords;
int field_length; 
int field_width; 
int common_field_length; 
int common_field_width; 
double* array_data;
double* dub_array_data;
double* dub_array;
double* check_array;
double* discrep_array; 
double* hpcG_array; 
double* discrep_lap_array;
double* hpcG_lap_array; 
double* array_pointer; 
int procs_rank;
int n_dims;
int* coords;
int* dub_coords;
int* dim_sizes;
int* periods;
	
int shift_source, shift_dest; 
MPI_Comm cart_comm;
MPI_Comm ColComm;
MPI_Comm RowComm;  
MPI_Status status; 
MPI_Request s_request;
MPI_Request r_request;

double* left_border;
double* right_border;
double* top_border; 
double* bottom_border;
double* SborderHor;
double* SborderVer;
double xl_pos;
double xr_pos;
double yt_pos;
double yb_pos;
double central_pos;
double mid_x_coords;
double mid_y_coords;

double discrep_value;
double hpcg_value;
double alpha_param;
double iter_param;

void create_coords();
double fill_border(const double x_coord, const double y_coord);
double fill_cell(const double x_coord, const double y_coord);
double fill_solut(const double x_coord, const double y_coord);
void exchange_function(double* source_array, int border_block);
double laplasian(double* cell_array, int x_pos, int y_pos, int border_block);
double dot_product(double* left_array, double* right_array);
	
void Initialization(int data_length, int data_width, int procs_rank, int p_size);
double* Iteration(const int iter_num); 
double Iter_accuracy();
double Total_accuracy();
void ClearAll();


//------------------------------------------------------------------------------------------------------

void create_coords() {
	for (int i = 0; i < global_length; ++i) {
		x_coords[i] = ((DEFAULT_END_X - DEFAULT_START_X) * (pow(1.0 + i / (global_length - 1.0), exponential_coeff) - 1.0) / (pow(2.0, exponential_coeff) - 1.0));
	}
	for (int i = 0; i < global_width; ++i) {
		y_coords[i] = ((DEFAULT_END_Y - DEFAULT_START_Y) * (pow(1.0 + i / (global_width - 1.0), exponential_coeff) - 1.0) / (pow(2.0, exponential_coeff) - 1.0));
	}
}

double fill_border(const double x_coord, const double y_coord) {
	return (1 + sin(x_coord * y_coord));
}

double fill_cell(const double x_coord, const double y_coord) {
	return (x_coord * x_coord + y_coord * y_coord)*sin(x_coord * y_coord);
}

double fill_solut(const double x_coord, const double y_coord) {
	return (1 + sin(x_coord * y_coord));
}

void exchange_function(double* source_array, int border_block) {
	int rank_source = 0;
	int	rank_dest = 0;
	// right border to left border.
	MPI_Cart_shift(cart_comm, 0, -1, &rank_source, &rank_dest);
	for (int i = 0; i < field_width; ++i) {
		SborderVer[i + border_block * field_width] = source_array[i * field_length];
	}
	MPI_Sendrecv(&SborderVer[border_block * field_width], field_width, MPI_DOUBLE, rank_dest, 2,
		&left_border[border_block * field_width], field_width, MPI_DOUBLE, rank_source, 2, cart_comm, &status);
	// left border to right border.
	MPI_Cart_shift(cart_comm, 0, 1, &rank_source, &rank_dest);
	for (int i = 0; i < field_width; ++i) {
		SborderVer[i + border_block * field_width] = source_array[(i + 1) * field_length - 1];
	}
	MPI_Sendrecv(&SborderVer[border_block * field_width], field_width, MPI_DOUBLE, rank_dest, 2,
		&right_border[border_block * field_width], field_width, MPI_DOUBLE, rank_source, 2, cart_comm, &status);	
	// top border to low border
	MPI_Cart_shift(cart_comm, 1, -1, &rank_source, &rank_dest);
	for (int i = 0; i < field_length; ++i) {
		SborderHor[i + border_block * field_length] = source_array[i];
	}
	MPI_Sendrecv(&SborderHor[border_block * field_length], field_length, MPI_DOUBLE, rank_dest, 1,
		&bottom_border[border_block * field_length], field_length, MPI_DOUBLE, rank_source, 1, cart_comm, &status);

	// low border to top border.
	MPI_Cart_shift(cart_comm, 1, 1, &rank_source, &rank_dest);
	for (int i = 0; i < field_length; ++i) {
		SborderHor[i + border_block * field_length] = source_array[field_length * (field_width - 1) + i];
	}
	MPI_Sendrecv(&SborderHor[border_block * field_length], field_length, MPI_DOUBLE, rank_dest, 1,
		&top_border[border_block * field_length], field_length, MPI_DOUBLE, rank_source, 1, cart_comm, &status);

}

//Laplasian
double laplasian(double* cell_array, int x_pos, int y_pos, int border_block) {
	double result = 0.0;
	if (x_pos == 0) {
		xl_pos = right_border[y_pos + field_width * border_block];
	}
	else {
		xl_pos = cell_array[y_pos * field_length + x_pos - 1];
	}
	if (x_pos == field_length - 1) {
		xr_pos = left_border[y_pos + field_width * border_block];
	}
	else {
		xr_pos = cell_array[y_pos * field_length + x_pos + 1];
	}
	if (y_pos == 0) {
		yt_pos = top_border[x_pos + field_length * border_block];
	}
	else {
		yt_pos = cell_array[(y_pos - 1) * field_length + x_pos];
	}
	if (y_pos == field_width - 1) {
		yb_pos = bottom_border[x_pos + field_length * border_block];
	}
	else {
		yb_pos = cell_array[(y_pos + 1) * field_length + x_pos];
	}
	central_pos  = cell_array[y_pos * field_length + x_pos];
	
	xl_step = (x_coords[x_bias + x_pos] - x_coords[x_bias + x_pos - 1]);
	xr_step = (x_coords[x_bias + x_pos + 1] - x_coords[x_bias + x_pos]);
	yt_step = (y_coords[y_bias + y_pos] - y_coords[y_bias + y_pos - 1]);
	yb_step = (y_coords[y_bias + y_pos + 1] - y_coords[y_bias + y_pos]);
	mid_x_coords = (xl_step + xr_step) / 2.0;
	mid_y_coords = (yt_step + yb_step) / 2.0;
	result = result + ((central_pos - xl_pos) / xl_step - (xr_pos - central_pos) / xr_step) / mid_x_coords;
	result = result + ((central_pos - yb_pos) / yt_step - (yt_pos - central_pos) / yb_step) / mid_y_coords;
	return result;
}

double dot_product(double* left_array, double* right_array) {
	double local_dot_product = 0;
	double global_dot_product = 0;
	for (int i = 0; i < field_width; ++i) {
		for (int j = 0; j < field_length; ++j) {
			local_dot_product += left_array[i * field_length + j] * right_array[i * field_length + j];
		}
	}
	MPI_Allreduce(&local_dot_product, &global_dot_product, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
	return global_dot_product;
}

void Initialization(int data_length, int data_width, int p_rank, int p_size) {
	n_dims = 2;
	global_length = data_length;
	global_width = data_width;
	dim_sizes = (int *)malloc(n_dims*sizeof(int));
	periods = (int *)malloc(n_dims*sizeof(int));
	coords = (int *)malloc(n_dims*sizeof(int));
	dub_coords = (int *)malloc(n_dims*sizeof(int));

	procs_rank = p_rank;
	int x_cart_size; x_cart_size = 1;
	int y_cart_size; y_cart_size = 0;
	while (x_cart_size < (trunc(sqrt(p_size)))) {
		x_cart_size *= 2;
	}
	y_cart_size = p_size / x_cart_size;

	dim_sizes[0] = x_cart_size;
	dim_sizes[1] = y_cart_size;
	periods[0] = periods[1] = 0;
	int reorder;
	reorder = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes, periods, reorder, &cart_comm);

	field_length = (data_length - data_length % x_cart_size) / x_cart_size;
	field_width = (data_width - data_width % y_cart_size) / y_cart_size;
	common_field_length = field_length;
	common_field_width = field_width;
	MPI_Cart_coords(cart_comm, procs_rank, n_dims, coords);
	if (coords[0] < data_length % x_cart_size) field_length++;
	if (coords[1] < data_width % y_cart_size) field_width++;

	int Subdims[2];

  	Subdims[0] = 0;  // Dimensionality fixing
  	Subdims[1] = 1;  // The presence of the given dimension in the subgrid
  	MPI_Cart_sub(cart_comm, Subdims, &RowComm);
  
	Subdims[0] = 1;
  	Subdims[1] = 0;
 	MPI_Cart_sub(cart_comm, Subdims, &ColComm);

	array_data = (double *)malloc(field_length * field_width*sizeof(double));
	dub_array_data = (double *)malloc(field_length * field_width*sizeof(double));
	dub_array = (double *)malloc(field_length * field_width*sizeof(double));
	check_array = (double *)malloc(field_length * field_width*sizeof(double));
	discrep_array = (double *)malloc(field_length * field_width*sizeof(double));
	hpcG_array = (double *)malloc(field_length * field_width*sizeof(double));
	discrep_lap_array = (double *)malloc(field_length * field_width*sizeof(double));
	hpcG_lap_array = (double *)malloc(field_length * field_width*sizeof(double));
	array_pointer = NULL;
	exponential_coeff = DEFAULT_COEFF;
	x_coords = (double *)malloc(global_length*sizeof(double));
	y_coords = (double *)malloc(global_length*sizeof(double));
	top_border = (double *)malloc(field_length * 3*sizeof(double));
	bottom_border = (double *)malloc(field_length * 3*sizeof(double));
	left_border = (double *)malloc(field_width * 3*sizeof(double));
	right_border = (double *)malloc(field_width * 3*sizeof(double));
	SborderHor = (double *)malloc(field_length * 3*sizeof(double));
	SborderVer = (double *)malloc(field_width * 3*sizeof(double));
	x_coord = 0;
	y_coord = 0;
	
	create_coords();
	x_bias = common_field_length * coords[0] + coords[0] * (coords[0] < (global_length % dim_sizes[0])) + (global_length % dim_sizes[0]) * (coords[0] >= (global_length % dim_sizes[0]));
	y_bias = common_field_width * coords[1] + coords[1] * (coords[1] < (global_width % dim_sizes[1])) + (global_width % dim_sizes[1]) * (coords[1] >= (global_width % dim_sizes[1]));
	//printf("%d \n",  x_bias);
	for (int i = 0; i < field_width; ++i) {
		for (int j = 0; j < field_length; ++j) {
			x_coord = x_coords[x_bias + j];
			y_coord = y_coords[y_bias + i]; 

			if ((x_coord == DEFAULT_START_X) || (x_coord == DEFAULT_END_X) || (y_coord == DEFAULT_START_Y) || (y_coord == DEFAULT_END_Y)) {
				array_data[field_length * i + j] = fill_border(x_coord, y_coord);
				//dub_array_data[field_length * i + j] = fill_border(x_coord, y_coord);
				//printf("%f ",array_data[field_length * i + j]);
			}
			else {
				array_data[field_length * i + j] = DEFAULT_VALUE_P;
			}
			check_array[field_length * i + j] = fill_solut(x_coord, y_coord);

			discrep_array[field_width * j + i] = 0;
			hpcG_array[field_width * j + i] = 0;
			discrep_lap_array[field_width * j + i] = 0;
			hpcG_lap_array[field_width * j + i] = 0;
		}
	}
}

void ResultCollection (double* pCMatrix, double* pCblock, int Size,int BlockSize) {
  double * pResultRow;// = new double [Size*BlockSize];
  pResultRow = (double *)malloc(Size * BlockSize*sizeof(double));

  for (int i=0; i<BlockSize; i++) {
    MPI_Gather( &pCblock[i*BlockSize], BlockSize, MPI_DOUBLE, 
      &pResultRow[i*Size], BlockSize, MPI_DOUBLE, 0, RowComm);
  }

  if (coords[1] == 0) {
    MPI_Gather(pResultRow, BlockSize*Size, MPI_DOUBLE, pCMatrix, 
      BlockSize*Size, MPI_DOUBLE, 0, ColComm);
  }
  free(pResultRow);
}

double* Iteration(const int iter_num) {
	double* arr;
	// Step 1 - calculate discrepancy array values and hpcG array values for future iterations.
	exchange_function(array_data, 0);
#pragma omp parallel for schedule(static)
	for (int i = 0 + 1 * (coords[1] == 0); i < field_width - 1 * (coords[1] == dim_sizes[1] - 1); ++i) {
		for (int j = 0 + 1 * (coords[0] == 0); j < field_length - 1 * (coords[0] == dim_sizes[0] - 1); ++j) {
			discrep_array[i * field_length + j] = laplasian(array_data, j, i, 0) -
				fill_cell(x_coords[x_bias + j], y_coords[y_bias + i]);
			if (iter_num == 0) {
				hpcG_array[i * field_length + j] = discrep_array[i * field_length + j];
				//printf("%f\n",discrep_array[i * field_length + j]);
			}
		}
	}

	// Step 2 - create discrep_lap_array and hpcG_lap_array ( last one only in case iter_num == 0 ), we'll need them anyway.
	exchange_function(discrep_array, 1);
#pragma omp parallel for schedule(static)
	for (int i = 0 + 1 * (coords[1] == 0); i < field_width - 1 * (coords[1] == dim_sizes[1] - 1); ++i) {
		for (int j = 0 + 1 * (coords[0] == 0); j < field_length - 1 * (coords[0] == dim_sizes[0] - 1); ++j) {
			discrep_lap_array[i * field_length + j] = laplasian(discrep_array, j, i, 1);
			if (iter_num == 0) {
				hpcG_lap_array[i * field_length + j] = discrep_lap_array[i * field_length + j];
			}
		}
	}
	// Step 3 - calculate new hpcG array values. Use dub_array for storage. This is needed only for HPCG iterations ( iter_num > 0 ).
	if (iter_num != 0) {
		alpha_param = dot_product(discrep_lap_array, hpcG_array) / dot_product(hpcG_lap_array, hpcG_array);
		if (coords[1] == 0) { 
#pragma omp parallel for schedule(static)
			for (int i = 0; i < field_length; ++i) dub_array[i] = 0; 
		}
		if (coords[1] == dim_sizes[1] - 1) { 
#pragma omp parallel for schedule(static)
			for (int i = 0; i < field_length; ++i) dub_array[field_length *  (field_width - 1) + i] = 0; 
		}
		if (coords[0] == 0) { 
#pragma omp parallel for schedule(static)
			for (int i = 0; i < field_width; ++i) dub_array[field_length * i] = 0; 
		}
		if (coords[0] == dim_sizes[0] - 1) { 
#pragma omp parallel for schedule(static)
			for (int i = 0; i < field_width; ++i) dub_array[field_length * i + (field_length - 1)] = 0; 
		}
#pragma omp parallel for schedule(static)
		for (int i = 0 + 1 * (coords[1] == 0); i < field_width - 1 * (coords[1] == dim_sizes[1] - 1); ++i) {
			for (int j = 0 + 1 * (coords[0] == 0); j < field_length - 1 * (coords[0] == dim_sizes[0] - 1); ++j) {
				dub_array[i * field_length + j] = discrep_array[i * field_length + j] - alpha_param * hpcG_array[i * field_length + j];
			}
		}
		array_pointer = hpcG_array;
		hpcG_array = dub_array;
		dub_array = array_pointer;
	}

	// Step 4 - calculate current hpcG_lap_array, if iteration isn't zero.
	if (iter_num != 0) {
		exchange_function(hpcG_array, 2);
#pragma omp parallel for schedule(static)
		for (int i = 0 + 1 * (coords[1] == 0); i < field_width - 1 * (coords[1] == dim_sizes[1] - 1); ++i) {
			for (int j = 0 + 1 * (coords[0] == 0); j < field_length - 1 * (coords[0] == dim_sizes[0] - 1); ++j) {
				hpcG_lap_array[i * field_length + j] = laplasian(hpcG_array, j, i, 2);
			}
		}
	}

	// Step 5 - calculate actual values for new iteration. Use dub_array for storage.
	if (iter_num != 0) {
		iter_param = dot_product(discrep_array, hpcG_array) / dot_product(hpcG_lap_array, hpcG_array);
	}
	else {
		iter_param = dot_product(discrep_array, discrep_array) / dot_product(discrep_lap_array, discrep_array);
	}
	for (int i = 0 + 1 * (coords[1] == 0); i < field_width - 1 * (coords[1] == dim_sizes[1] - 1); ++i) {
		for (int j = 0 + 1 * (coords[0] == 0); j < field_length - 1 * (coords[0] == dim_sizes[0] - 1); ++j) {
			if (iter_num != 0) {
				hpcG_lap_array[i * field_length + j] = laplasian(hpcG_array, j, i, 2);
				dub_array[i * field_length + j] = array_data[i * field_length + j] - iter_param * hpcG_array[i * field_length + j];
				//printf("%f \n",array_data[field_length + j]);
				}
			else {
				discrep_lap_array[i * field_length + j] = laplasian(discrep_array, j, i, 1);
				dub_array[i * field_length + j] = array_data[i * field_length + j]- iter_param * discrep_array[i * field_length + j];
				//printf("%f \n",dub_array[field_length + j]);

				}

		}
	}
	arr = dub_array;	

	array_pointer = array_data;
	array_data = dub_array;
	dub_array = array_pointer;
	return arr;
}

double Iter_accuracy() {
	for (int i = 0; i < field_width; ++i) {
		for (int j = 0; j < field_length; ++j) {
			dub_array[i * field_length + j] = array_data[i * field_length + j] - dub_array[i * field_length + j];
		}
	}
	return sqrt(dot_product(dub_array, dub_array) / (global_length * global_width));
}

double Total_accuracy() {
	for (int i = 0; i < field_width; ++i) {
		for (int j = 0; j < field_length; ++j) {
			check_array[field_length * i + j] = check_array[field_length * i + j] - array_data[field_length * i + j];
		}
	}
	return sqrt(dot_product(check_array, check_array) / (global_length * global_width));
}

void ClearAll(){
	if (dim_sizes != NULL) { free(dim_sizes);
	}
	if (periods != NULL) { free(periods);
	}
	if (dub_coords != NULL) { free(dub_coords);
	}
	if (coords != NULL) { free(coords);
	}
	if (array_data != NULL) { free(array_data);
	}
	if (dub_array != NULL) { free(dub_array);
	}
	if (discrep_array != NULL) { free(discrep_array);
	}
	if (hpcG_array != NULL) {free(hpcG_array);
	}
	if (discrep_lap_array != NULL) { free(discrep_lap_array);
	}
	if (hpcG_lap_array != NULL) { free(hpcG_lap_array);
	}
	if (check_array != NULL) { free(check_array);
	}
	if (top_border != NULL) { free(top_border);
	}
	if (bottom_border != NULL) { free(bottom_border);
	}
	if (left_border != NULL) { free(left_border);
	}
	if (right_border != NULL) { free(right_border);
	}
	if (SborderHor != NULL) { free(SborderHor);
	}
	if (SborderVer != NULL) { free(SborderVer);
	}
}