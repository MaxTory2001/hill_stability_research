/**
 * Hill Stability
 * 
 * Here I use the IAS15 integrator to check stability of orbits
 * in heirarchical triples for different inclinations and semi-major axes
 * The integrator automatically adjusts the timestep so that 
 * even very high eccentricity encounters are resolved with high 
 * accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include<stdbool.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
void print_time(struct reb_simulation* r);
void generate_sims(struct reb_simulation* r, float inc, float a_rh);
double tmax = 2.*M_PI*1.0e2;     // 100 years translated into code units
bool stable;

int main(int argc, char* argv[]){
    system("rm -v stuff_alternate.txt");        // delete previous output file

    // Initial conditions
    for (float inc = 0.0; inc <= 180.0; inc++){
        for (float a_over_rh = 0.3; a_over_rh < 1.1; a_over_rh += 0.01){
        
            struct reb_simulation* r = reb_create_simulation();
            // Setup constants
            r->dt             = 1e-2;     // initial timestep
            r->integrator        = REB_INTEGRATOR_WHFAST;
            r->heartbeat        = heartbeat;

            // Adding star
            struct reb_particle star = {0}; 
            star.m  = 1;
            star.hash = reb_hash("star");
            reb_add(r, star); 

            generate_sims(r, inc * M_PI / 180, a_over_rh);

            reb_free_simulation(r);
        }
    }
}

void generate_sims(struct reb_simulation* r, float inc, float a_rh){
    double r_hill = pow((1.01e-6)/3., 1.0/3.0);
    double a = a_rh * r_hill;
    double v_c = sqrt(1.0e-6/a);
    // the planet (m1)
    struct reb_particle m1 = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0e-6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    m1.hash = reb_hash("m1");

    reb_add(r, m1);

    // The orbiting satellite (m2)
    struct reb_particle m2 = {0};
    // trying to calculate manuallyy - is this different?
    m2.m = 1e-8;
    m2.x = m1.x + a;
    m2.vx = m1.vx;
    m2.vy = m1.vy + v_c * cos(inc);
    m2.vz = m1.vz + v_c * sin(inc);
    m2.hash = reb_hash("m2"); 
    
    reb_add(r, m2); 

    reb_move_to_com(r);

    printf("\ninc = %.4f, a = %.4f\n", inc, a_rh);

    stable = true;

    reb_integrate(r, tmax);

    if (stable) print_time(r);
}

void heartbeat(struct reb_simulation* r){
    float dist = sqrt(pow(r->particles[1].x - r->particles[2].x, 2) + pow(r->particles[1].y - r->particles[2].y, 2) + pow(r->particles[1].z - r->particles[2].z, 2));
    float r_hill = pow((r->particles[1].m + r->particles[2].m)/3, 1.0/3.0);  // not sure about this one chief

    if(reb_output_check(r, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(r, tmax);
    }
    if(dist >= r_hill * 3.0 && stable){            // outputs to a file
        print_time(r);
        stable = false;
    }
}

void print_time(struct reb_simulation* r){
    FILE* fp = fopen("stuff_alternate.txt", "a");
    fprintf(fp, "%.4f\n", r->t);
    fclose(fp);
}
