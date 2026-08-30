/*
 * lorenz_chaos.c
 * --------------
 * Simulates the Lorenz system (the origin of "the butterfly effect") using
 * a 4th-order Runge-Kutta integrator, and numerically estimates its largest
 * Lyapunov exponent -- the number that quantifies HOW chaotic the system is.
 *
 * The Lorenz system:
 *   dx/dt = sigma * (y - x)
 *   dy/dt = x * (rho - z) - y
 *   dz/dt = x * y - beta * z
 *
 * With the classic parameters (sigma=10, rho=28, beta=8/3) this system is
 * chaotic: bounded, deterministic, but with trajectories that diverge
 * exponentially from arbitrarily close starting points. That divergence
 * rate is the Lyapunov exponent lambda. If lambda > 0, the system is chaotic;
 * the value of lambda tells you the "doubling time" of the divergence
 * (roughly, prediction horizon ~ 1/lambda).
 *
 * This program does two things:
 *   1. Integrates ONE trajectory and dumps it to lorenz_trajectory.csv
 *      (for plotting the classic butterfly-shaped attractor).
 *   2. Integrates TWO trajectories starting a tiny distance epsilon apart,
 *      periodically rescaling the separation back to epsilon (the standard
 *      Benettin algorithm), and accumulates the average exponential growth
 *      rate of their separation -> the largest Lyapunov exponent.
 *
 * Build:  gcc -O2 -o lorenz_chaos lorenz_chaos.c -lm
 * Run:    ./lorenz_chaos
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Lorenz system parameters (the classic chaotic regime) */
#define SIGMA 10.0
#define RHO   28.0
#define BETA  (8.0/3.0)

typedef struct { double x, y, z; } State;

static State lorenz_deriv(State s) {
    State d;
    d.x = SIGMA * (s.y - s.x);
    d.y = s.x * (RHO - s.z) - s.y;
    d.z = s.x * s.y - BETA * s.z;
    return d;
}

static State state_add(State a, State b, double scale) {
    State r = { a.x + b.x * scale, a.y + b.y * scale, a.z + b.z * scale };
    return r;
}

/* One step of classic 4th-order Runge-Kutta */
static State rk4_step(State s, double dt) {
    State k1 = lorenz_deriv(s);
    State k2 = lorenz_deriv(state_add(s, k1, dt / 2.0));
    State k3 = lorenz_deriv(state_add(s, k2, dt / 2.0));
    State k4 = lorenz_deriv(state_add(s, k3, dt));

    State sum = { k1.x + 2*k2.x + 2*k3.x + k4.x,
                  k1.y + 2*k2.y + 2*k3.y + k4.y,
                  k1.z + 2*k2.z + 2*k3.z + k4.z };
    return state_add(s, sum, dt / 6.0);
}

static double dist(State a, State b) {
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

int main(void) {
    const double dt = 0.005;
    const int transient_steps = 2000;   /* let the trajectory settle onto the attractor first */
    const int trajectory_steps = 8000;  /* steps saved for plotting */

    /* ---------- Part 1: single trajectory for plotting the attractor ---------- */
    State s = { 1.0, 1.0, 1.0 };
    for (int i = 0; i < transient_steps; i++) s = rk4_step(s, dt);

    FILE *traj = fopen("lorenz_trajectory.csv", "w");
    fprintf(traj, "t,x,y,z\n");
    for (int i = 0; i < trajectory_steps; i++) {
        fprintf(traj, "%.5f,%.6f,%.6f,%.6f\n", i * dt, s.x, s.y, s.z);
        s = rk4_step(s, dt);
    }
    fclose(traj);
    fprintf(stderr, "Wrote lorenz_trajectory.csv (%d points)\n", trajectory_steps);

    /* ---------- Part 2: Benettin algorithm for the largest Lyapunov exponent ---------- */
    /* Two trajectories start epsilon apart; every `renorm_every` steps we measure
     * how much they diverged, accumulate log(growth), then rescale the perturbed
     * trajectory back to distance epsilon along the same direction. Averaging
     * log(growth)/dt over many renormalizations converges to the Lyapunov exponent. */
    const double epsilon = 1e-8;
    const int renorm_every = 10;
    const int n_renorms = 4000;

    State ref = { 1.0, 1.0, 1.0 };
    for (int i = 0; i < transient_steps; i++) ref = rk4_step(ref, dt);

    State pert = ref;
    pert.x += epsilon; /* perturb slightly in x */

    double sum_log_growth = 0.0;
    FILE *lyap_log = fopen("lyapunov_convergence.csv", "w");
    fprintf(lyap_log, "renorm_step,running_lambda\n");

    for (int r = 0; r < n_renorms; r++) {
        for (int i = 0; i < renorm_every; i++) {
            ref  = rk4_step(ref, dt);
            pert = rk4_step(pert, dt);
        }
        double d = dist(ref, pert);
        sum_log_growth += log(d / epsilon);

        /* rescale perturbed trajectory back to distance epsilon from ref,
         * preserving direction of divergence */
        double scale = epsilon / d;
        pert.x = ref.x + (pert.x - ref.x) * scale;
        pert.y = ref.y + (pert.y - ref.y) * scale;
        pert.z = ref.z + (pert.z - ref.z) * scale;

        double elapsed_time = (r + 1) * renorm_every * dt;
        double running_lambda = sum_log_growth / elapsed_time;
        fprintf(lyap_log, "%d,%.6f\n", r, running_lambda);
    }
    fclose(lyap_log);

    double total_time = n_renorms * renorm_every * dt;
    double lambda = sum_log_growth / total_time;

    fprintf(stderr, "\n=== Largest Lyapunov exponent estimate ===\n");
    fprintf(stderr, "lambda ~= %.4f  (per unit time)\n", lambda);
    fprintf(stderr, "Known reference value for classic Lorenz params: ~0.90-0.91\n");
    fprintf(stderr, "Prediction horizon (1/lambda): ~%.3f time units\n", 1.0 / lambda);
    fprintf(stderr, "lambda > 0  =>  system is chaotic: nearby trajectories diverge exponentially.\n");

    return 0;
}
