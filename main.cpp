#include <boost/thread.hpp>
#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/assert.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <random>

using namespace std;

// Function Declarations
void init_simulation(vector<vector<double>> &atom_pos, int atoms_per_side, double atom_spacing, double L);
void internal_E(vector<vector<double>> &pos, double cutt_off, int N, double L, double *E);
void E_to_file(fstream &file, double E);
void calc_P(vector<vector<double>> &pos, double cutt_off, int N, double T, double rho, double L, double *P);
void run_simulation(vector<vector<double>> atom_pos, double rc, double L, double temperature, int N_mcs_eq,
                    int N_mcs, double *E_cumulative, double *P_cumulative, int *num_accpted);


int main() {
    int N = 200; //number of atoms
    double rho = 0.9; //density
    double rc; //potential cut off distance
    double E_cumulative; // Cumulative E
    double P_cumulative; // Cumulative P
    double acc_ratio; // for calculating acceptance ratios
    double temperature = 2.0;
    int N_mcs = 500; // Number of data collection Monte Carlo Steps
    int N_mcs_eq = 500; // Number of Equlibration Monte Carlo Steps

    // Find the smallest cube such that we have at minimum the dersired number of atoms
    int atoms_per_side = 2;
    while (atoms_per_side * atoms_per_side * atoms_per_side < N) atoms_per_side++;
    N = atoms_per_side * atoms_per_side * atoms_per_side; // reset N

    vector<vector<double>> atom_pos(N, vector<double>(3)); // Vector of atom positions


    // Set up cubic lattice and L such that we have to correct density and number of atoms
    double L = cbrt((double)(N)/ rho);
    double atom_spacing = L / (double) atoms_per_side;
    rc = L/ 2.0;
    init_simulation(atom_pos, atoms_per_side, atom_spacing, L); // Builds the cubic lattice


    cout << "----------------------- LJ MC SIMULATION ------------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "Running simulation with the following parameters:" << endl;
    cout << "N: " << N << endl;
    cout << "L: " << L << endl;
    cout << "rho: " << rho << endl;
    cout << "cutt off: " << rc << endl;
    cout << "T: " << temperature << endl;
    cout << "Equlibration steps: " << (N_mcs_eq * N) << endl;
    cout << "Data Collection Steps: " << (N_mcs * N) << endl;
    cout << endl;
    cout << endl;
    cout << "----------------------- STARTING SIMULATION ----------------------" << endl;
    cout << endl;
    cout << endl;

    // Set up worker thread variables
    double E_cummulative_1 = 0.0;
    double E_cummulative_2 = 0.0;
    double E_cummulative_3 = 0.0;
    double E_cummulative_4 = 0.0;

    double P_cummulative_1 = 0.0;
    double P_cummulative_2 = 0.0;
    double P_cummulative_3 = 0.0;
    double P_cummulative_4 = 0.0;

    int num_acc_1 = 0;
    int num_acc_2 = 0;
    int num_acc_3 = 0;
    int num_acc_4 = 0;

    // Launch the worker Threads
    boost::timer t; // timer for the simulation
    boost::thread worker1(run_simulation, atom_pos, rc, L, temperature, N_mcs_eq, N_mcs/4, &E_cummulative_1, &P_cummulative_1, &num_acc_1);
    boost::thread worker2(run_simulation, atom_pos, rc, L, temperature, N_mcs_eq, N_mcs/4, &E_cummulative_2, &P_cummulative_2, &num_acc_2);
    boost::thread worker3(run_simulation, atom_pos, rc, L, temperature, N_mcs_eq, N_mcs/4, &E_cummulative_3, &P_cummulative_3, &num_acc_3);
    boost::thread worker4(run_simulation, atom_pos, rc, L, temperature, N_mcs_eq, N_mcs/4, &E_cummulative_4, &P_cummulative_4, &num_acc_4);

    // Wait for worker Threads to finish
    worker1.join();
    worker2.join();
    worker3.join();
    worker4.join();

    // Work out all the averages and print to console
    E_cumulative = (E_cummulative_1/(double)(N_mcs)) + (E_cummulative_2/(double)(N_mcs)) + (E_cummulative_3/(double)(N_mcs)) + (E_cummulative_4/(double)(N_mcs));
    P_cumulative = (P_cummulative_1/(double)(N_mcs)) + (P_cummulative_2/(double)(N_mcs)) + (P_cummulative_3/(double)(N_mcs)) + (P_cummulative_4/(double)(N_mcs));
    acc_ratio = (double)(num_acc_1 + num_acc_2 + num_acc_3 + num_acc_4) / (double)(N_mcs * N);

    // Calculate the pressure and Energy corrections
    double r_correction = 1.0 / (rc * rc * rc);
    double E_correction = (8.0 / 3.0) * M_PI * rho * ((r_correction * r_correction * r_correction) / 3.0 - r_correction);
    double P_correction = (16.0 / 3.0) * M_PI * rho * rho * ((2.0 / 3.0) * r_correction * r_correction * r_correction - r_correction);


    cout << "----------------------- SIMULATION RESULTS -----------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "Time taken: " << t.elapsed()/4.0 << " Seconds" << endl;
    cout << "Acceptance ratio: " << acc_ratio << endl;
    cout << "Average Total Energy: " << E_cumulative + (N * E_correction) << endl;
    cout << "Average Energy per paticle: " << E_cumulative/N + E_correction << endl;
    cout << "Pressure Correction: " << P_correction << endl;
    cout << "Average Total Pressure: " << P_cumulative + P_correction << endl;
    cout << endl;
    cout << endl;
    cout << "----------------------- END OF RESULTS --------------------------" << endl;

    return 0;
}

// Function to calculate the internal energy of the system using a LJ potentail
void internal_E(vector<vector<double>> &pos, double cutt_off, int N, double L, double *E) {
    double x;
    double y;
    double z;
    double r;

    *E = 0.0;

    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {

            // Find Difference in atom positions
            x = pos[i][0] - pos[j][0];
            y = pos[i][1] - pos[j][1];
            z = pos[i][2] - pos[j][2];

            // Apply minimum image convention
            if (x > L / 2.0) x -= L;
            else if (x < -L / 2.0) x +=L;
            if (y > L / 2.0) y -= L;
            else if (y < -L / 2.0) y +=L;
            if (z > L / 2.0) z -= L;
            else if (z < -L / 2.0) z +=L;

            // Assert the minimum image convention
            BOOST_ASSERT_MSG(fabs(x) <= L / 2.0, "Check the minimum image convention X > L/2");
            BOOST_ASSERT_MSG(fabs(y) <= L / 2.0, "Check the minimum image convention Y > L/2");
            BOOST_ASSERT_MSG(fabs(z) <= L / 2.0, "Check the minimum image convention Z > L/2");

            // find r^2 distance
            r = x * x + y * y + z * z;

            if (r > cutt_off * cutt_off) { // check r^2 !> rc^2
                continue; // if it is just increment the counter and go back to the top
            }

            // compute (1/r)^6 did this to avoid using Pow as this is quicker and has less side effects
            r = 1.0 / (r * r * r);

            // Add the contribution to the energy to the total energy
            *E += 4.0 * ((r*r) - r);
        }
    }
}

// Simple function to dump energys to a file used to check for equlibated system
void E_to_file(fstream &file, double E) {
    if (file.is_open()){
        file << E << endl;
    }
    else {
        cout << "FILE ISNT OPEN" << endl;
    }
}

// Function to calculate the pressure using the virial as per frenkel et al.
void calc_P(vector<vector<double>> &pos, double cutt_off, int N, double T, double rho, double L, double *P) {
    double x;
    double y;
    double z;
    double r;
    double vir = 0.0;
    double v = L * L * L;

    // Compute the virial sum
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {

            // Find difference in position
            x = pos[i][0] - pos[j][0];
            y = pos[i][1] - pos[j][1];
            z = pos[i][2] - pos[j][2];

            // Apply minimum image convention
            if (x > L / 2.0) x -= L;
            else if (x < -L / 2.0) x +=L;
            if (y > L / 2.0) y -= L;
            else if (y < -L / 2.0) y +=L;
            if (z > L / 2.0) z -= L;
            else if (z < -L / 2.0) z +=L;

            // Assert the minimum image convention
            BOOST_ASSERT_MSG(fabs(x) <= L / 2.0, "Check the minimum image convention X > L/2");
            BOOST_ASSERT_MSG(fabs(y) <= L / 2.0, "Check the minimum image convention Y > L/2");
            BOOST_ASSERT_MSG(fabs(z) <= L / 2.0, "Check the minimum image convention Z > L/2");

            // compute r^2
            r = x * x + y * y + z * z;

            // Check we are within the cutoff
            if (r > cutt_off * cutt_off) {
                continue;
            }

            // Compute (1/r)^6
            r = 1.0 / (r * r * r);

            // Sum the contribution
            vir += 48.0 * (r*r - 0.5*r);
        }
    }

    vir /= 3.0;

    // compute the actual pressure using the above virial sum
    *P = (rho * T) + (vir / v);
}

// Function to initalise a cubic lattice
void init_simulation(vector<vector<double>> &atom_pos, int atoms_per_side, double atom_spacing, double L) {

    int atom_count = 0;
    for (int i = 0; i < atoms_per_side; ++i) {
        for (int j = 0; j < atoms_per_side; ++j) {
            for (int k = 0; k < atoms_per_side; ++k) {
                atom_pos[atom_count][0] = fmod((double)k * atom_spacing, atom_spacing * atoms_per_side); // x position
                atom_pos[atom_count][1] = fmod((double)j * atom_spacing, atom_spacing * atoms_per_side); // y position
                atom_pos[atom_count][2] = fmod((double)i * atom_spacing, atom_spacing * atoms_per_side); // z position

                // Assertions to check the atoms are actually in the box useful for debugging but turned off in release compile
                BOOST_ASSERT_MSG(atom_pos[atom_count][0] <= L && atom_pos[atom_count][0] >= 0.0, "atom initialised out of box in x direction");
                BOOST_ASSERT_MSG(atom_pos[atom_count][1] <= L && atom_pos[atom_count][1] >= 0.0, "atom initialised out of box in y direction");
                BOOST_ASSERT_MSG(atom_pos[atom_count][2] <= L && atom_pos[atom_count][2] >= 0.0, "atom initialised out of box in z direction");

                atom_count++;
            }
        }

    }
}

// Function to run an MC simulation
void run_simulation(vector<vector<double>> atom_pos, double rc, double L, double temperature, int N_mcs_eq, int N_mcs, double *E_cumulative, double *P_cumulative, int *num_accpted) {
    double E_new;
    double P_new;
    double E_curr;
    double P_curr;
    int acc = 0;
    int choice;
    int N = (int)atom_pos.size();
    double lambda = 0.1; // larges possible move
    vector<double> rand_vec(3);
    vector<double> old_pos(3);
    random_device rd; // used to seed the random number generators such that each thread has its own indipendent stream of RN's

    // Sort out the PRNG'S
    boost::random::mt19937 gen(rd());
    boost::uniform_int<> rand_int(0, N-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<>> rand_choice(gen, rand_int);
    boost::uniform_on_sphere<double> unit_vec(3);
    boost::variate_generator<boost::mt19937&, boost::uniform_on_sphere<double >> random_on_sphere(gen, unit_vec);
    boost::uniform_01<double> uni_01;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<double>> met_rand(gen, uni_01);

    // Calculate the initial E and P
    internal_E(atom_pos, rc, N, L, &E_curr);
    calc_P(atom_pos, rc, N, temperature, (double)N/(L*L*L), L, &P_curr);

    // Equlibrate the system
    for (int l = 0; l < N_mcs_eq; ++l) { // each monte carlo step
        for (int i = 0; i < N; ++i) { // each atom
            choice = rand_choice(); // choose random atom
            rand_vec = random_on_sphere(); //random unit vector

            // save old position
            old_pos = atom_pos[choice];

            // Move atom
            atom_pos[choice][0] += lambda * rand_vec[0];
            atom_pos[choice][1] += lambda * rand_vec[1];
            atom_pos[choice][2] += lambda * rand_vec[2];

            // make sure the atoms stay in the box
            atom_pos[choice][0] += L;
            atom_pos[choice][1] += L;
            atom_pos[choice][2] += L;

            atom_pos[choice][0] = fmod(atom_pos[choice][0],L);
            atom_pos[choice][1] = fmod(atom_pos[choice][1],L);
            atom_pos[choice][2] = fmod(atom_pos[choice][2],L);


            // ASSERT the new positions useful for debugging but turned of during release build
            BOOST_ASSERT_MSG(atom_pos[choice][0] <= L && atom_pos[choice][0] >= 0.0, "atom moved out of box in x direction");
            BOOST_ASSERT_MSG(atom_pos[choice][1] <= L && atom_pos[choice][1] >= 0.0, "atom moved out of box in x direction");
            BOOST_ASSERT_MSG(atom_pos[choice][2] <= L && atom_pos[choice][2] >= 0.0, "atom moved out of box in x direction");

            // Calc new E and P
            internal_E(atom_pos, rc, N, L, &E_new);
            calc_P(atom_pos, rc, N, temperature, (double)N/(L*L*L), L, &P_new);

            // Test Metropolis condition
            if (met_rand() < exp(-(E_new - E_curr)/ temperature)) {
                // Accept Move and set E and P
                E_curr = E_new;
                P_curr = P_new;
            }
            else {
                // Reject and move atoms back
                atom_pos[choice] = old_pos;
            }
        }
    }

    // Record Data

    for (int m = 0; m < N_mcs; ++m) { //each monte carlo step
        for (int i = 0; i < N; ++i) { // each atom
            choice = rand_choice(); // choose random atom
            rand_vec = random_on_sphere(); //random unit

            // save old position
            old_pos = atom_pos[choice];

            // Move atom
            atom_pos[choice][0] += lambda * rand_vec[0];
            atom_pos[choice][1] += lambda * rand_vec[1];
            atom_pos[choice][2] += lambda * rand_vec[2];

            // make sure the atoms stay in the box
            atom_pos[choice][0] += L;
            atom_pos[choice][1] += L;
            atom_pos[choice][2] += L;

            atom_pos[choice][0] = fmod(atom_pos[choice][0],L);
            atom_pos[choice][1] = fmod(atom_pos[choice][1],L);
            atom_pos[choice][2] = fmod(atom_pos[choice][2],L);


            // ASSERT the new positions useful for debugging but turned of during release build
            BOOST_ASSERT_MSG(atom_pos[choice][0] <= L && atom_pos[choice][0] >= 0.0, "atom moved out of box in x direction");
            BOOST_ASSERT_MSG(atom_pos[choice][1] <= L && atom_pos[choice][1] >= 0.0, "atom moved out of box in x direction");
            BOOST_ASSERT_MSG(atom_pos[choice][2] <= L && atom_pos[choice][2] >= 0.0, "atom moved out of box in x direction");

            // Calc new E and P
            internal_E(atom_pos, rc, N, L, &E_new);
            calc_P(atom_pos, rc, N, temperature, (double)N/(L*L*L), L, &P_new);

            // Test Metropolis condition
            if (met_rand() < exp(-(E_new - E_curr)/ temperature)) {
                // Accept Move and E and P
                E_curr = E_new;
                P_curr = P_new;
                acc++;
            }
            else {
                // Reject and move atoms back
                atom_pos[choice] = old_pos;
            }
        }
        // After Each MCS collect data for the averaging
        *E_cumulative += E_curr;
        *P_cumulative += P_curr;
    }
    *num_accpted = acc;
}