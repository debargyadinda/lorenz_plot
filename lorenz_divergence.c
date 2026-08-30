

#include <stdio.h>
#include <math.h>

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
    const int transient_steps = 2000;
    const int steps = 6000; /* ~30 time units: long enough for full dramatic divergence */
    const double epsilon = 1e-8;

    State a = { 1.0, 1.0, 1.0 };
    for (int i = 0; i < transient_steps; i++) a = rk4_step(a, dt);

    State b = a;
    b.x += epsilon;

    FILE *f = fopen("lorenz_divergence.csv", "w");
    fprintf(f, "t,x_a,y_a,z_a,x_b,y_b,z_b,separation\n");
    for (int i = 0; i < steps; i++) {
        double t = i * dt;
        double sep = dist(a, b);
        fprintf(f, "%.5f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.10f\n",
                t, a.x, a.y, a.z, b.x, b.y, b.z, sep);
        a = rk4_step(a, dt);
        b = rk4_step(b, dt);
    }
    fclose(f);
    fprintf(stderr, "Wrote lorenz_divergence.csv (%d points, starting separation %.0e)\n",
            steps, epsilon);
    return 0;
}
