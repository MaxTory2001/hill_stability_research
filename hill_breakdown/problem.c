/**
 * Hill Stability Breakdown
 * 
 * Here I use the IAS15 integrator to check stability of orbits
 * in heirarchical triples for different inclinations and semi-major axes,
 * then extend this to less hierarchical regimes to determine the breakdown
 * point of the hill stability criterion
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
void print_time(struct reb_simulation* r);
void generate_sims(struct reb_simulation* r);
void global_orbit_output(struct reb_simulation* r, struct reb_particle p);
void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p);
void loop_sma();

double tmax = 2.*M_PI*1.0e2;     // 100 years translated into code units
bool stable;
int n = 0; float inc; float a_over_rh; float mass_ratio = 0.0;

// For plotting info about a specific orbit
int target_inc = 50;
float target_a = 0.4;

// hill radius
float r_hill;

int main(int argc, char* argv[]){
    system("rm -v time*.txt");        // delete previous output file
    system("rm -v diagnostics.txt");        // delete previous output file
    system("rm -v detailed_global.txt");
    system("rm -v detailed_relative.txt"); 

    // Initial conditions
    if (argc > 1){
        // Mass ratio between planet and outer star
        mass_ratio = pow(10, strtof(argv[1], NULL));

        r_hill = pow((mass_ratio)/3., 1.0/3.0);

        if (argc == 3) {
            // Only check one given inclination - for computation time
            inc = strtof(argv[2], NULL);
            loop_sma();
        
        } else {
            for (inc = 0.0; inc <= 180.0; inc++){
                loop_sma();
            }
        } 
    } else {
        printf("Must provide mass ratio as command line argument!\n");
        return 1;
    }

    
}

void loop_sma(){
    /*
     *
    */
    for (a_over_rh = 0.2; a_over_rh < 1.5; a_over_rh += 0.01){
        
        struct reb_simulation* r = reb_create_simulation();
        // Setup constants
        r->dt             = M_PI*1.0e-3;     // initial timestep
        r->integrator        = REB_INTEGRATOR_IAS15;
        r->heartbeat        = heartbeat;

        // Adding star
        struct reb_particle planet = {0}; 
        planet.m  = mass_ratio;
        planet.hash = reb_hash("m1");
        reb_add(r, planet); 

        generate_sims(r);

        // 1 in ~100 orbits have their diagnostics saved. 97 is prime so we get a good combination of inclinations and a's
        if (n % 97 == 0){
            FILE* fp = fopen("diagnostics.txt", "a");
            fprintf(fp, "\n");
            fclose(fp);
        }

        n++;

        reb_free_simulation(r);
    }
}

void generate_sims(struct reb_simulation* r){
    // the satelite (m2)
    struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0e-2 * r->particles[0].m, a_over_rh * r_hill, 0.0, inc * M_PI / 180, 0.0, 0.0, 0.0);
    m2.hash = reb_hash("m2");

    reb_add(r, m2);

    // The star
    struct reb_particle star = reb_tools_orbit_to_particle(r->G, r->particles[1],  1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    star.hash = reb_hash("star"); 

    reb_add(r, star); 

    reb_move_to_com(r);

    printf("\ninc = %.4f, a = %.4f, n = %d\n", inc, a_over_rh, n);

    reb_integrate(r, tmax);
    // orbit was not disrupted, output time
    print_time(r);
}

void heartbeat(struct reb_simulation* r){
    float dist = reb_particle_distance(&r->particles[0], &r->particles[1]);
    /*
    // detailed information about the orbit specified
    if (n == target_inc * 81 + round(100 * (target_a-0.3))){
        global_orbit_output(r, r->particles[1]);
        relative_orbit_output(r, r->particles[1], r->particles[0]);
    }

    // diagnostics about 1 in every 97 orbits
    if(reb_output_check(r, 5*M_PI)){        // outputs to the screen
        if (n % 97 == 0){
            FILE* fp = fopen("diagnostics.txt", "a");
            struct reb_vec3d ang_mom = reb_tools_angular_momentum(r);
            fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, reb_tools_energy(r), ang_mom.x, ang_mom.y, ang_mom.z);
            fclose(fp);
        }
    }*/
    
    if(reb_output_check(r, 200*M_PI)){        // outputs to the screen
        reb_output_timing(r, tmax);
    }

    // printing when each orbit is disrupted
    if(dist >= r_hill * 3.0 && stable){            
        r->status = 1;
    }
}

void print_time(struct reb_simulation* r){
    char fname[20];

    sprintf(fname, "time_%.2e.txt", mass_ratio);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.7f\n", r->t);
    fclose(fp);
}

void global_orbit_output(struct reb_simulation* r, struct reb_particle p) {
    // prints detailed info about an orbit - orbital elements and more
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, p, r->particles[2]);
    FILE* fp = fopen("detailed_global.txt", "a");
    fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, p.x, p.y, p.z, p.vx, p.vy, p.vz, o.a, o.e, o.inc, o.Omega);
    fclose(fp);
}

void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p) {
    /*
    * p: primary particle - centre of the orbit
    * s: secondary particle
    */
    // prints detailed info about an orbit - orbital elements and more
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, s, p);
    FILE* fp = fopen("detailed_relative.txt", "a");
    fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, (s.x-p.x) / r_hill, (s.y-p.y) / r_hill, (s.z-p.z) / r_hill, s.vx-p.vx, s.vy-p.vy, s.vz-p.vz, o.a, o.e, o.inc, o.omega);
    fclose(fp);
}
