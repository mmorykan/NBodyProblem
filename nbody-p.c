/**
 * Runs a simulation of the n-body problem in 3D.
 * 
 * To compile the program:
 *   gcc -Wall -O3 -march=native -fopenmp nbody-p.c matrix.c util.c -o nbody-p -lpthread -lm
 * 
 * To run the program:
 *   ./nbody-p time-step total-time outputs-per-body input.npy output.npy [opt: num-threads]
 * where:
 *   - time-step is the amount of time between steps (Δt, in seconds)
 *   - total-time is the total amount of time to simulate (in seconds)
 *   - outputs-per-body is the number of positions to output per body
 *   - input.npy is the file describing the initial state of the system (below)
 *   - output.npy is the output of the program (see below)
 *   - last argument is an optional number of threads (a reasonable default is
 *     chosen if not provided)
 * 
 * input.npy has a n-by-7 matrix with one row per body and the columns:
 *   - mass (in kg)
 *   - initial x, y, z position (in m)
 *   - initial x, y, z velocity (in m/s)
 * 
 * output.npy is generated and has a (T/Δt)-by-(3n) matrix with each row
 * containing the x, y, and z positions of each of the n bodies after a given
 * timestep.
 * 
 * See the PDF for implementation details and other requirements.
 * 
 * AUTHORS: Jonah Beers and Mark Morykan
 */

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

#include "matrix.h"
#include "util.h"
#include "helpers.h"

/**
 * Store all positions for each body for each time step into an output matrix
 */
void compute_n_bodies(nbody* info) {
    size_t output_steps = info->output_steps, output_rows = info->output_rows, num_threads = info->num_threads,
            num_bodies = info->num_bodies, num_steps = info->num_steps;
    double* curr_pos = info->curr_pos, *next_pos = info->next_pos;
    dbl_rstrct_pntr masses = info->masses;
    dbl_rstrct_pntr curr_vel = info->curr_vel;
    dbl_rstrct_pntr output = info->output;
    double time_step = info->time_step;

    #pragma omp parallel num_threads(num_threads)
    {
        // Compute every position for every body for every time step
        for (size_t time = 0; time < num_steps; time++) {
            size_t index = time/output_steps*output_rows;  // place in output matrix
            // Each thread gets a different chunk of bodies to compute
            #pragma omp for
            for (size_t i = 0; i < num_bodies; i++) {
                size_t x_i = i * 3, y_i = x_i + 1, z_i = y_i + 1;  // places in position and velocity arrays
                double acc[3] = {0, 0, 0};  // reset acceleration array for each body
                double x = curr_pos[x_i], y = curr_pos[y_i], z = curr_pos[z_i];
                for (size_t j = 0; j < num_bodies; j++) {
                    if (i == j) continue;  // can't compute attraction between a body and itself
                    store_acceleration(acc, curr_pos, masses[j], x_i, y_i, z_i, j);  // store values in the acceleration array
                }

                acc[0] *= G, acc[1] *= G, acc[2] *= G;  // multiply accelerations by Gravitational constant
                
                // Store new positions into next_pos by computing new velocities
                double velocity_x = curr_vel[x_i] += acc[0] * time_step;
                double velocity_y = curr_vel[y_i] += acc[1] * time_step;
                double velocity_z = curr_vel[z_i] += acc[2] * time_step;
                
                next_pos[x_i] = x + velocity_x * time_step;
                next_pos[y_i] = y + velocity_y * time_step;
                next_pos[z_i] = z + velocity_z * time_step; 

                store_output(output, index, x, y, z, x_i, y_i, z_i, time, output_steps);  // store positions in the output matrix
            }
            // make sure only one thread swaps the current position with the next position
            #pragma omp single
            {
                swap(&curr_pos, &next_pos); // make the next position array the current position array for next time step's calculations
            }
        }
    }
}


int main(int argc, const char* argv[]) {
    // Variables available now:
    //   time_step    number of seconds between each time point
    //   total_time   total number of seconds in the simulation
    //   num_steps    number of time steps to simulate (more useful than total_time)
    //   num_outputs  number of times the position will be output for all bodies
    //   output_steps number of steps between each output of the position
    //   num_threads  number of threads to use
    //   input        n-by-7 Matrix of input data
    //   n            number of bodies to simulate

    // Parse arguments
    if (argc != 6 && argc != 7) { fprintf(stderr, "usage: %s time-step total-time outputs-per-body input.npy output.npy [num-threads]\n", argv[0]); return 1; }
    double time_step = atof(argv[1]), total_time = atof(argv[2]);
    if (time_step <= 0 || total_time <= 0 || time_step > total_time) { fprintf(stderr, "time-step and total-time must be positive with total-time > time-step\n"); return 1; }
    size_t num_outputs = atoi(argv[3]);
    if (num_outputs <= 0) { fprintf(stderr, "outputs-per-body must be positive\n"); return 1; }
    size_t num_threads = argc == 7 ? atoi(argv[6]) : get_num_cores_affinity()/2; 
    if (num_threads <= 0) { fprintf(stderr, "num-threads must be positive\n"); return 1; }
    Matrix* input = matrix_from_npy_path(argv[4]);
    if (input == NULL) { perror("error reading input"); return 1; }
    if (input->cols != 7) { fprintf(stderr, "input.npy must have 7 columns\n"); return 1; }
    size_t n = input->rows;
    if (n == 0) { fprintf(stderr, "input.npy must have at least 1 row\n"); return 1; }
    if (num_threads > n) { num_threads = n; }
    size_t num_steps = (size_t)(total_time / time_step + 0.5);
    if (num_steps < num_outputs) { num_outputs = 1; }
    size_t output_steps = num_steps/num_outputs;
    num_outputs = (num_steps+output_steps-1)/output_steps;

    // Start the clock
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    size_t output_rows = 3*n;
    Matrix* output = matrix_zeros(num_outputs, output_rows);  // output matrix

    // Copy all input positions, masses, and velocities into arrays to maximize cache hits
    size_t position_space = output_rows*sizeof(double);
    double* masses   = (double*)malloc(n*sizeof(double));
    double* curr_pos = (double*)malloc(position_space);
    double* next_pos = (double*)malloc(position_space);
    double* curr_vel = (double*)malloc(position_space);
    store_matrix_info(input->data, input->cols, n, masses, curr_pos, curr_vel);

    // Set up struct with all body info 
    nbody* info = (nbody*)malloc(sizeof(nbody));
    info->output = output->data;
    info->output_steps = output_steps;
    info->output_rows = output_rows;
    info->num_bodies = n;
    info->masses = masses;
    info->curr_pos = curr_pos;
    info->next_pos = next_pos;
    info->curr_vel = curr_vel;
    info->time_step = time_step;
    info->num_steps = num_steps;
    info->num_threads = num_threads;

    compute_n_bodies(info);  // main portion of the program

    // Get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    // Save results
    matrix_to_npy_path(argv[5], output);

    // Cleanup
    matrix_free(input);
    matrix_free(output);
    free(masses);
    free(curr_pos);
    free(next_pos);
    free(curr_vel);
    free(info);
    
    return 0;
}
