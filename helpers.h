#include <stdlib.h>

// Gravitational Constant in N m^2 / kg^2 or m^3 / kg / s^2
#define G 6.6743015e-11

// Softening factor to reduce divide-by-near-zero effects
#define SOFTENING 1e-9

// Double* restrict
#define dbl_rstrct_pntr double* restrict

// Used to send info to compute_n_bodies function
typedef struct {
    size_t num_bodies, output_steps, output_rows, num_steps, num_threads;
    double *masses, *curr_pos, *next_pos, *curr_vel, *input_data, *output;
    double time_step;
} nbody;

/**
 * Store masses, velocities, and positions into separate arrays from the input matrix
 */
void store_matrix_info(double* input_data, size_t input_cols, size_t num_bodies, double* masses, double* curr_pos, double* curr_vel) {
    for (size_t i = 0; i < num_bodies; i++) {
        size_t curr_pos_row = 3*i, input_row = i*input_cols;

        masses[i] = input_data[input_row];

        curr_pos[curr_pos_row] = input_data[input_row+1];
        curr_pos[curr_pos_row+1] = input_data[input_row+2];
        curr_pos[curr_pos_row+2] = input_data[input_row+3];

        curr_vel[curr_pos_row] = input_data[input_row+4];
        curr_vel[curr_pos_row+1] = input_data[input_row+5];
        curr_vel[curr_pos_row+2] = input_data[input_row+6];
    }
}

/**
 * Store result positions in output matrix
 */ 
void store_output(double* output, size_t index, double x, double y, double z, size_t x_i, size_t y_i, size_t z_i, size_t time, size_t output_steps) {
    if (time % output_steps == 0) {
        output[index + x_i] = x;
        output[index + y_i] = y;
        output[index + z_i] = z;
    }
}

/**
 * Make the current position array the next position array
 */
void swap(double** curr_pos, double** next_pos) {
    double* temp = *curr_pos;
    *curr_pos = *next_pos;
    *next_pos = temp;
}

/**
 * Stores a body's acceleration
 */
inline static void store_acceleration(double* acc, double* curr_pos, double mass, size_t x_i, size_t y_i, size_t z_i, size_t j) {
    size_t x_j = j * 3, y_j = x_j + 1, z_j = y_j + 1;
    double x_diff = curr_pos[x_j] - curr_pos[x_i], y_diff = curr_pos[y_j] - curr_pos[y_i], z_diff = curr_pos[z_j] - curr_pos[z_i];
    double dist = sqrt((x_diff * x_diff) + (y_diff * y_diff) + (z_diff * z_diff)) + SOFTENING;
    double temp = mass / (dist * dist * dist);

    acc[0] += x_diff * temp;  
    acc[1] += y_diff * temp;
    acc[2] += z_diff * temp;
}
