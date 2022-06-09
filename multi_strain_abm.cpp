// Extracts Python objects from an ABMlux state file, converts these objects to arrays of C data
// types and executes a basic version of the ABMlux multi-strain model

#include "Python.h"
#include <iostream>
#include <stdlib.h>
#include <ctime>

// g++ -o baseline_strains .\baseline_strains.cpp -L .\libs\*

// Psuedo-random number generator seed
#define PRNG_SEED 0

// Clock
#define SIMULATION_LENGTH_TICKS 57600

// State file
char state_file[] = "states/state_strains_25.abm";

// Constants
#define TICKS_PER_WEEK 1008
#define MAX_DISEASE_PROFILE_LENGTH 8
#define NUMBER_OF_HEALTH_STATES 9
#define NUMBER_OF_STRAINS 2
#define MAX_LOCATIONS_FOR_ACTIVITY 10
#define NUMBER_OF_ACTIVITIES 14
#define HOUSE_ACTIVITY 0
#define PUBLIC_TRANSPORT_ACTIVITY 6

// Health states
#define SUSCEPTIBLE 0
#define EXPOSED 1
#define ASYMPTOMATIC 2
#define PREINFECTED 3
#define INFECTED 4
#define HOSPITALIZING 5
#define VENTILATING 6
#define RECOVERED 7
#define DEAD 8

// Interventions
#define LOCKDOWN_START 3024
#define LOCKDOWN_END 7056

// Iterator
#define INF_INDEX_SIZE 24
#define TRANS_INDEX_SIZE 24
double incubation_period_multiplier[NUMBER_OF_STRAINS] = {1.0, 1.0};
double infectious_period_multiplier[INF_INDEX_SIZE][NUMBER_OF_STRAINS]
    = {{1.0, 0.2}, {1.0, 0.4}, {1.0, 0.6}, {1.0, 0.8}, {1.0, 1.0}, {1.0, 1.2},
       {1.0, 1.4}, {1.0, 1.6}, {1.0, 1.8}, {1.0, 2.0}, {1.0, 2.2}, {1.0, 2.4},
       {1.0, 2.6}, {1.0, 2.8}, {1.0, 3.0}, {1.0, 3.2}, {1.0, 3.4}, {1.0, 3.6},
       {1.0, 3.8}, {1.0, 4.0}, {1.0, 4.2}, {1.0, 4.4}, {1.0, 4.6}, {1.0, 4.8}};
double transmission_probability_multiplier[TRANS_INDEX_SIZE][NUMBER_OF_STRAINS]
    = {{1.0, 0.2}, {1.0, 0.4}, {1.0, 0.6}, {1.0, 0.8}, {1.0, 1.0}, {1.0, 1.2},
       {1.0, 1.4}, {1.0, 1.6}, {1.0, 1.8}, {1.0, 2.0}, {1.0, 2.2}, {1.0, 2.4},
       {1.0, 2.6}, {1.0, 2.8}, {1.0, 3.0}, {1.0, 3.2}, {1.0, 3.4}, {1.0, 3.6},
       {1.0, 3.8}, {1.0, 4.0}, {1.0, 4.2}, {1.0, 4.4}, {1.0, 4.6}, {1.0, 4.8}};

int main(){

    // Extracts Python simulator object, simulates baseline and lockdown scenarios and calculates
    // the difference in total cumulative cases, iterated over a range of parameter multipliers

    // ################################# EXTRACT SIMULATOR #########################################

    Py_Initialize();

    // Extract simulator
    PyObject* integer_model_filename
        = PyUnicode_FromString((char*)"abm_sim_extractor");
    PyObject* integer_model
        = PyImport_Import(integer_model_filename);
    PyObject* extract_sim
        = PyObject_GetAttrString(integer_model, (char*)"extract_sim");
    PyObject* args
        = PyTuple_Pack(1, PyUnicode_FromString((char*)state_file));
    PyObject* sim
        = PyObject_CallObject(extract_sim, args);

    // Extract agents
    PyObject* build_integer_agents
        = PyObject_GetAttrString(integer_model, (char*)"build_integer_agents");
    PyObject* integer_agents
        = PyObject_CallObject(build_integer_agents, PyTuple_Pack(1, sim));
    int number_of_agents = PyLong_AsLong(PyTuple_GetItem(integer_agents, 0));
    auto resident = new int[number_of_agents]();
    for (int n = 0; n < number_of_agents; ++n){
        resident[n] = PyLong_AsLong(PyList_GetItem(PyTuple_GetItem(integer_agents, 1), n));
    }
    int rescale_factor = PyLong_AsLong(PyTuple_GetItem(integer_agents, 2));

    // Extract locations
    PyObject* build_integer_locations
        = PyObject_GetAttrString(integer_model, (char*)"build_integer_locations");
    PyObject* integer_locations
        = PyObject_CallObject(build_integer_locations, PyTuple_Pack(1, sim));
    int number_of_locations = PyLong_AsLong(PyTuple_GetItem(integer_locations, 0));
    int num_hospitals = PyLong_AsLong(PyTuple_GetItem(integer_locations, 1));
    auto hospitals = new int[num_hospitals]();
    for (int n = 0; n < num_hospitals; ++n){
        hospitals[n] = PyLong_AsLong(PyList_GetItem(PyTuple_GetItem(integer_locations, 2), n));
    }
    int num_cemeteries = PyLong_AsLong(PyTuple_GetItem(integer_locations, 3));
    auto cemeteries = new int[num_cemeteries]();
    for (int n = 0; n < num_cemeteries; ++n){
        cemeteries[n] = PyLong_AsLong(PyList_GetItem(PyTuple_GetItem(integer_locations, 4), n));
    }
    int num_pt_locations = PyLong_AsLong(PyTuple_GetItem(integer_locations, 5));
    auto pt_locations = new int[num_pt_locations]();
    for (int n = 0; n < num_pt_locations; ++n){
        pt_locations[n] = PyLong_AsLong(PyList_GetItem(PyTuple_GetItem(integer_locations, 6), n));
    }
    auto pt_availability = new int[TICKS_PER_WEEK]();
    for (int n = 0; n < TICKS_PER_WEEK; ++n){
        pt_availability[n]
            = PyLong_AsLong(PyList_GetItem(PyTuple_GetItem(integer_locations, 7), n));
    }
    auto initial_location = new int[number_of_agents]();
    for (int n = 0; n < number_of_agents; ++n){
        initial_location[n]
            = PyLong_AsLong(PyList_GetItem(PyTuple_GetItem(integer_locations, 8), n));
    }

    // Extract activities
    PyObject* build_integer_activities
        = PyObject_GetAttrString(integer_model, (char*)"build_integer_activities");
    PyObject* integer_activities
        = PyObject_CallObject(build_integer_activities, PyTuple_Pack(1, sim));
    auto weekly_routine = new int[number_of_agents][TICKS_PER_WEEK]();
    auto num_locations_for_activity = new int[number_of_agents][NUMBER_OF_ACTIVITIES]();
    auto locations_for_activity
        = new int[number_of_agents][NUMBER_OF_ACTIVITIES][MAX_LOCATIONS_FOR_ACTIVITY]();
    for (int n = 0; n < number_of_agents; ++n){
        for (int t = 0; t < TICKS_PER_WEEK; ++t){
            weekly_routine[n][t] = PyLong_AsLong(PyList_GetItem(PyList_GetItem(PyTuple_GetItem(
                integer_activities, 0), n), t));
        }
        for (int a = 0; a < NUMBER_OF_ACTIVITIES; ++a){
            num_locations_for_activity[n][a] = PyLong_AsLong(PyList_GetItem(PyList_GetItem(
                PyTuple_GetItem(integer_activities, 1), n), a));
            for (int l = 0; l < MAX_LOCATIONS_FOR_ACTIVITY; ++l){
                locations_for_activity[n][a][l] = PyLong_AsLong(PyList_GetItem(PyList_GetItem(
                    PyList_GetItem(PyTuple_GetItem(integer_activities, 2), n), a), l));
            }
        }
    }

    // Extract disease
    PyObject* build_integer_disease
        = PyObject_GetAttrString(integer_model, (char*)"build_integer_disease");
    PyObject* integer_disease
        = PyObject_CallObject(build_integer_disease, PyTuple_Pack(1, sim));
    auto transmission_probabilities = new double[NUMBER_OF_STRAINS][2]();
    auto disease_durations_dict
        = new int[number_of_agents][NUMBER_OF_STRAINS][MAX_DISEASE_PROFILE_LENGTH]();
    auto num_initial_cases_by_strain = new int[NUMBER_OF_STRAINS]();
    for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
        transmission_probabilities[i][0] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(
            PyTuple_GetItem(integer_disease, 0), i), 0));
        transmission_probabilities[i][1] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(
            PyTuple_GetItem(integer_disease, 0), i), 1));
        num_initial_cases_by_strain[i] = PyLong_AsLong(PyList_GetItem(
            PyTuple_GetItem(integer_disease, 4), i));
    }
    auto transmission_multiplier = new double[number_of_locations]();
    for (int m = 0; m < number_of_locations; ++m){
        transmission_multiplier[m]
            = PyFloat_AsDouble(PyList_GetItem(PyTuple_GetItem(integer_disease, 1), m));
    }
    auto disease_profile_dict
        = new int[number_of_agents][NUMBER_OF_STRAINS][MAX_DISEASE_PROFILE_LENGTH]();
    for (int n = 0; n < number_of_agents; ++n){
        for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
            for (int d = 0; d < MAX_DISEASE_PROFILE_LENGTH; ++d){
                disease_profile_dict[n][i][d] = PyLong_AsLong(PyList_GetItem(PyList_GetItem(
                    PyList_GetItem(PyTuple_GetItem(integer_disease, 2), n), i), d));
                disease_durations_dict[n][i][d] = PyLong_AsLong(PyList_GetItem(PyList_GetItem(
                    PyList_GetItem(PyTuple_GetItem(integer_disease, 3), n), i), d));
            }
        }
    }

    // ################################# SIMULATE ##################################################

    srand(PRNG_SEED);

    // Initialize agent health
    printf("Generating initial health state...\n");
    auto initial_strain = new int[number_of_agents]();
    for (int n = 0; n < number_of_agents; ++n){initial_strain[n] = -1;}
    for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
        int selected_cases;
        selected_cases = 0;
        while (selected_cases < num_initial_cases_by_strain[i]){
            int ind;
            ind = (int)(number_of_agents * ((double)((rand() * RAND_MAX) + rand()) /
                (double)((RAND_MAX * RAND_MAX) + 1)));
            if (initial_strain[ind] == -1){
                if (resident[ind] == 1){initial_strain[ind] = i; selected_cases += 1;}
            }
        }
    }

    // Initialize final output
    int resident_total_cumulative_cases[INF_INDEX_SIZE][TRANS_INDEX_SIZE][2];
    for (int inf_multiplier_index = 0;
         inf_multiplier_index < INF_INDEX_SIZE;
         ++inf_multiplier_index){
        for (int trans_multiplier_index = 0;
             trans_multiplier_index < TRANS_INDEX_SIZE;
             ++trans_multiplier_index){
            resident_total_cumulative_cases[inf_multiplier_index][trans_multiplier_index][0] = 0;
            resident_total_cumulative_cases[inf_multiplier_index][trans_multiplier_index][1] = 0;
        }
    }

    // For a given location and strain, the array counts[][] will record the number of agents in
    // that location infected with that strain, distinguishing between those who are in health
    // states ASYMPTOMATIC or PREINFECTED from those who are in health states INFECTED,
    // HOSPITALIZING or VENTILATING, since the transmission probabilities may differ. For a given
    // location m, the matrices mat_a[m] and mat_b[m] are used to calculate the probability that an
    // agent, located in m, becomes infected with a particular strain. The vector weight[m] is also
    // used when calculating the above probability.

    auto counts = new int[number_of_locations][NUMBER_OF_STRAINS][2]();
    auto mat_a = new double[number_of_locations][NUMBER_OF_STRAINS][2]();
    auto mat_b = new double[number_of_locations][NUMBER_OF_STRAINS][2]();
    auto weight = new double[number_of_locations][NUMBER_OF_STRAINS]();

    // Iterate over multipliers
    for (int inf_multiplier_index = 0;
         inf_multiplier_index < INF_INDEX_SIZE;
         ++inf_multiplier_index){
        for (int trans_multiplier_index = 0;
             trans_multiplier_index < TRANS_INDEX_SIZE;
             ++trans_multiplier_index){

            // Apply multipliers
            for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
                transmission_probabilities[i][0] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(
                    PyTuple_GetItem(integer_disease, 0), i), 0));
                transmission_probabilities[i][1] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(
                    PyTuple_GetItem(integer_disease, 0), i), 1));
                transmission_probabilities[i][0] *=
                    transmission_probability_multiplier[trans_multiplier_index][i];
                transmission_probabilities[i][1] *=
                    transmission_probability_multiplier[trans_multiplier_index][i];
                for (int n = 0; n < number_of_agents; ++n){
                    for (int d = 0; d < MAX_DISEASE_PROFILE_LENGTH; ++d){
                        disease_durations_dict[n][i][d] =
                            PyLong_AsLong(PyList_GetItem(PyList_GetItem(PyList_GetItem(
                                PyTuple_GetItem(integer_disease, 3), n), i), d));
                        int health_state;
                        health_state = disease_profile_dict[n][i][d];
                        if (health_state == EXPOSED){
                            disease_durations_dict[n][i][d] = (int)(disease_durations_dict[n][i][d]
                                                                 * incubation_period_multiplier[i]);
                        }
                        if (health_state == ASYMPTOMATIC ||
                            health_state == PREINFECTED ||
                            health_state == INFECTED ||
                            health_state == HOSPITALIZING ||
                            health_state == VENTILATING){
                            disease_durations_dict[n][i][d] = (int)(disease_durations_dict[n][i][d]
                                           * infectious_period_multiplier[inf_multiplier_index][i]);
                        }
                    }
                }
            }

            // Iterate over baseline and lockdown scenarios for these mutlipliers
            for (int scenario = 0; scenario < 2; ++scenario){

                // Reset agent current state
                int *current_location = new int[number_of_agents];
                for (int n = 0; n < number_of_agents; ++n){
                    current_location[n] = initial_location[n];
                }
                int *current_health = new int[number_of_agents];
                for (int n = 0; n < number_of_agents; ++n){current_health[n] = 0;}
                int *current_strain = new int[number_of_agents];
                for (int n = 0; n < number_of_agents; ++n){current_strain[n] = -1;}
                int *disease_profile_index_dict = new int[number_of_agents];
                for (int n = 0; n < number_of_agents; ++n){disease_profile_index_dict[n] = 0;}

                // Reset PRNG (note rand() has 'resolution' approx 0.00003)
                srand(PRNG_SEED);

                // Precalculate the matrices mat_a and mat_b
                for (int m = 0; m < number_of_locations; ++m){
                    for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
                        mat_a[m][i][0] = transmission_multiplier[m]
                                         * transmission_probabilities[i][0];
                        mat_a[m][i][1] = transmission_multiplier[m]
                                         * transmission_probabilities[i][1];
                        mat_b[m][i][0] = 1 - mat_a[m][i][0];
                        mat_b[m][i][1] = 1 - mat_a[m][i][1];
                    }
                }

                // Reset health state counts
                int *resident_health_state_counts = new int[NUMBER_OF_HEALTH_STATES];
                int *resident_strain_counts = new int[NUMBER_OF_STRAINS];
                int *resident_cumulative_strain_counts = new int[NUMBER_OF_STRAINS];
                for (int r = 0; r < NUMBER_OF_HEALTH_STATES; ++r){
                    resident_health_state_counts[r] = 0;
                }
                for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
                    resident_strain_counts[i] = 0;
                    resident_cumulative_strain_counts[i] = 0;
                }

                // Reset health state change times
                int *health_state_change_time = new int[number_of_agents];
                for (int n = 0; n < number_of_agents; ++n){
                    health_state_change_time[n] = 0;
                }

                // Reset health states
                for (int n = 0; n < number_of_agents; ++n){
                    current_strain[n] = initial_strain[n];
                    if (current_strain[n] != -1){
                        disease_profile_index_dict[n] = 2;
                        current_health[n] = disease_profile_dict[n][current_strain[n]][2];
                    } else {
                        current_health[n] = SUSCEPTIBLE;
                    }
                }

                // Perform initial counts
                for (int m = 0; m < number_of_locations; ++m){
                    for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
                        counts[m][i][0] = 0;
                        counts[m][i][1] = 0;
                    }
                }
                for (int n = 0; n < number_of_agents; ++n){
                    if (resident[n] == 1){resident_health_state_counts[current_health[n]] += 1;}
                    if (current_health[n] == ASYMPTOMATIC ||
                        current_health[n] == PREINFECTED){
                        counts[current_location[n]][current_strain[n]][0] += 1;
                    }
                    if (current_health[n] == INFECTED ||
                        current_health[n] == HOSPITALIZING ||
                        current_health[n] == VENTILATING){
                        counts[current_location[n]][current_strain[n]][1] += 1;
                    }
                }

                // Start clock
                std::clock_t start;
                start = std::clock();
                printf("Starting simulation scenario %d inf %.1f trans %.1f...\n",
                       scenario, infectious_period_multiplier[inf_multiplier_index][1],
                       transmission_probability_multiplier[trans_multiplier_index][1]);

                // Initialize iteration output
                FILE *fpt_health;
                char filename_health_counts[32];
                sprintf(filename_health_counts, "output/resident_health_counts_scenario_"
                                                "%d_inf_%.1f_trans_%.1f_prng_%d.csv",
                        scenario, infectious_period_multiplier[inf_multiplier_index][1],
                        transmission_probability_multiplier[trans_multiplier_index][1], PRNG_SEED);
                fpt_health = fopen(filename_health_counts, "w+");
                fprintf(fpt_health,"SUSCEPTIBLE, EXPOSED, ASYMPTOMATIC, PREINFECTED, INFECTED,"
                                   " HOSPITALIZING, VENTILATING, RECOVERED, DEAD\n");
                FILE *fpt_strain;
                char filename_strain_counts[32];
                sprintf(filename_strain_counts, "output/resident_strain_counts_scenario_"
                                                "%d_inf_%.1f_trans_%.1f_prng_%d.csv",
                        scenario, infectious_period_multiplier[inf_multiplier_index][1],
                        transmission_probability_multiplier[trans_multiplier_index][1], PRNG_SEED);
                fpt_strain = fopen(filename_strain_counts, "w+");
                FILE *fpt_cumula;
                char filename_cumulative_counts[32];
                sprintf(filename_cumulative_counts, "output/resident_cumulative_counts_scenario_"
                                                    "%d_inf_%.1f_trans_%.1f_prng_%d.csv",
                        scenario, infectious_period_multiplier[inf_multiplier_index][1],
                        transmission_probability_multiplier[trans_multiplier_index][1], PRNG_SEED);
                fpt_cumula = fopen(filename_cumulative_counts, "w+");

                // Simulate
                for (int t = 0; t < SIMULATION_LENGTH_TICKS; ++t){

                    // Print ticks
                    printf("\r%d / %d", t, SIMULATION_LENGTH_TICKS - 1);
                    fflush(stdout);

                    // Calculate time offsets
                    int t_now, t_next;
                    t_now = t % (TICKS_PER_WEEK);
                    t_next = (t+1) % (TICKS_PER_WEEK);

                    // Initialize update arrays
                    int *health_updates = new int[number_of_agents];
                    int *needs_health_update = new int[number_of_agents];
                    int *strain_updates = new int[number_of_agents];
                    int *location_updates = new int[number_of_agents];
                    int *needs_location_update = new int[number_of_agents];
                    for (int n = 0; n < number_of_agents; ++n){
                        health_updates[n] = 0;
                        needs_health_update[n] = 0;
                        strain_updates[n] = 0;
                        location_updates[n] = 0;
                        needs_location_update[n] = 0;
                    }

                    // Calculate transmission weights for each location in integer form
                    double *weight_sum = new double[number_of_locations];
                    long *p_long = new long[number_of_locations];
                    for (int m = 0; m < number_of_locations; ++m){
                        double p_double = 1.0;
                        weight_sum[m] = 0.0;
                        for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
                            p_double *= pow(mat_b[m][i][0], (double) counts[m][i][0]);
                            p_double *= pow(mat_b[m][i][1], (double) counts[m][i][1]);
                            weight[m][i] = (mat_a[m][i][0] * (double) counts[m][i][0])
                                           + (mat_a[m][i][1] * (double) counts[m][i][1]);
                            weight_sum[m] += weight[m][i];
                        }
                        p_long[m] = (long)((1.0 - p_double) * (double) ((RAND_MAX * RAND_MAX) + 1));
                    }

                    // Loop through the agents
                    for (int n = 0; n < number_of_agents; ++n){

                        // Get current health and location of agent
                        int health, location;
                        health = current_health[n];
                        location = current_location[n];

                        // Determine health transitions
                        if (health == SUSCEPTIBLE){
                            if (p_long[location] > 0){
                                int r_1, r_2;
                                r_1 = rand();
                                r_2 = rand();
                                if((long)((r_1 * RAND_MAX) + r_2) < p_long[location]){
                                    health_updates[n] = EXPOSED;
                                    double rnd;
                                    rnd = ((double) rand()) * weight_sum[location];
                                    int strain = 0;
                                    while (rnd >= ((double) (RAND_MAX + 1))
                                                  * weight[location][strain]){
                                        rnd -= ((double) (RAND_MAX + 1))
                                               * weight[location][strain];
                                        ++strain;
                                    }
                                    strain_updates[n] = strain;
                                    needs_health_update[n] = 1;
                                }
                            }
                        }
                        if (health == EXPOSED ||
                            health == ASYMPTOMATIC ||
                            health == PREINFECTED ||
                            health == INFECTED ||
                            health == HOSPITALIZING ||
                            health == VENTILATING){
                            if (t - health_state_change_time[n] > disease_durations_dict[n][
                                current_strain[n]][disease_profile_index_dict[n]]){
                                health_updates[n] = disease_profile_dict[n][
                                    current_strain[n]][disease_profile_index_dict[n] + 1];
                                needs_health_update[n] = 1;
                            }
                        }

                        // Determine location transitions
                        if (weekly_routine[n][t_next] != weekly_routine[n][t_now]){
                            if (health != HOSPITALIZING && health != VENTILATING && health != DEAD){
                                int random_int_loc, new_activity, new_location, num_locs;
                                random_int_loc = rand();
                                new_activity = weekly_routine[n][t_next];
                                if (new_activity == PUBLIC_TRANSPORT_ACTIVITY){
                                    num_locs = (int) fmin((float) num_pt_locations,
                                                          (float) pt_availability[t_next]);
                                    new_location = pt_locations[(int)(((float) random_int_loc /
                                                   (float) (RAND_MAX + 1)) * (float) num_locs)];
                                    location_updates[n] = new_location;
                                    needs_location_update[n] = 1;
                                } else {
                                    num_locs = num_locations_for_activity[n][new_activity];
                                    int index_rand = (int)(((float) random_int_loc /
                                                     (float) (RAND_MAX + 1)) * (float) num_locs);
                                    new_location
                                        = locations_for_activity[n][new_activity][index_rand];
                                    location_updates[n] = new_location;
                                    needs_location_update[n] = 1;
                                }
                                if (scenario == 1){
                                    if (t >= LOCKDOWN_START && t < LOCKDOWN_END){
                                        location_updates[n]
                                            = locations_for_activity[n][HOUSE_ACTIVITY][0];
                                    }
                                }
                            }
                        }
                    }

                    delete[] p_long;
                    delete[] weight_sum;

                    // Update agent states
                    for (int n = 0; n < number_of_agents; ++n){

                        // Apply health updates
                        if (needs_health_update[n] == 1){
                            current_health[n] = health_updates[n];
                            disease_profile_index_dict[n] += 1;
                            health_state_change_time[n] = t;
                            needs_health_update[n] = 0;

                            // Add or remove strain
                            if (health_updates[n] == EXPOSED){
                                current_strain[n] = strain_updates[n];
                                if (resident[n] == 1){
                                    resident_cumulative_strain_counts[current_strain[n]] += 1;
                                }
                            }
                            if (health_updates[n] == RECOVERED || health_updates[n] == DEAD){
                                current_strain[n] = -1;
                            }
                        }

                        // Apply location updates
                        if (needs_location_update[n] == 1){
                            current_location[n] = location_updates[n];
                            needs_location_update[n] = 0;
                        }
                    }

                    delete[] health_updates;
                    delete[] needs_health_update;
                    delete[] strain_updates;
                    delete[] location_updates;
                    delete[] needs_location_update;

                    // Hospitalization intervention
                    for (int n = 0; n < number_of_agents; ++n){
                        if (current_health[n] == HOSPITALIZING || current_health[n] == VENTILATING){
                            int in_hospital = 0;
                            for (int h = 0; h < num_hospitals; ++h){
                                if (current_location[n] == hospitals[h]){
                                    in_hospital = 1;
                                }
                            }
                            if (in_hospital == 0){
                                int random_int_hosp, hospital;
                                random_int_hosp = rand();
                                hospital = hospitals[(int)(((float) random_int_hosp /
                                    (float) (RAND_MAX + 1)) * (float) num_hospitals)];
                                current_location[n] = hospital;
                            }
                        }
                        if (current_health[n] == DEAD){
                            int in_cemetery = 0;
                            for (int h = 0; h < num_cemeteries; ++h){
                                if (current_location[n] == cemeteries[h]){
                                    in_cemetery = 1;
                                }
                            }
                            if (in_cemetery == 0){
                                int random_int_ceme, cemetery;
                                random_int_ceme = rand();
                                cemetery = cemeteries[(int)(((float) random_int_ceme /
                                    (float) (RAND_MAX + 1)) * (float) num_cemeteries)];
                                current_location[n] = cemetery;
                            }
                        }
                    }

                    // Recount
                    for (int m = 0; m < number_of_locations; ++m){
                        for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
                            counts[m][i][0] = 0;
                            counts[m][i][1] = 0;
                        }
                    }
                    for (int i = 0; i < NUMBER_OF_STRAINS; ++i){
                        resident_strain_counts[i] = 0;
                    }
                    for (int r = 0; r < NUMBER_OF_HEALTH_STATES; ++r){
                        resident_health_state_counts[r] = 0;
                    }
                    for (int n = 0; n < number_of_agents; ++n){
                        if (resident[n] == 1){
                            resident_health_state_counts[current_health[n]] += 1;
                            if (current_strain[n] != -1){
                                resident_strain_counts[current_strain[n]] += 1;
                            }
                        }
                        if (current_health[n] == ASYMPTOMATIC ||
                            current_health[n] == PREINFECTED){
                            counts[current_location[n]][current_strain[n]][0] += 1;
                        }
                        if (current_health[n] == INFECTED ||
                            current_health[n] == HOSPITALIZING ||
                            current_health[n] == VENTILATING){
                            counts[current_location[n]][current_strain[n]][1] += 1;
                        }
                    }

                    // Prinf iteration output to files
                    fprintf(fpt_health,"%d, %d, %d, %d, %d, %d, %d, %d, %d\n",
                            rescale_factor * resident_health_state_counts[0],
                            rescale_factor * resident_health_state_counts[1],
                            rescale_factor * resident_health_state_counts[2],
                            rescale_factor * resident_health_state_counts[3],
                            rescale_factor * resident_health_state_counts[4],
                            rescale_factor * resident_health_state_counts[5],
                            rescale_factor * resident_health_state_counts[6],
                            rescale_factor * resident_health_state_counts[7],
                            rescale_factor * resident_health_state_counts[8]);
                    fprintf(fpt_strain,"%d, %d\n", rescale_factor * resident_strain_counts[0],
                                                   rescale_factor * resident_strain_counts[1]);
                    fprintf(fpt_cumula,"%d, %d\n",
                            rescale_factor * resident_cumulative_strain_counts[0],
                            rescale_factor * resident_cumulative_strain_counts[1]);

                }

                // Close iteration output files
                fclose(fpt_health);
                fclose(fpt_strain);
                fclose(fpt_cumula);

                // Print simulation duration
                printf("\n");
                printf("Simulation finished successfully in %f seconds\n",
                       (std::clock() - start) / (double) CLOCKS_PER_SEC);

                // Record total cumulative cases for this scenario for these multipliers
                resident_total_cumulative_cases[inf_multiplier_index][
                    trans_multiplier_index][scenario]= rescale_factor *
                    (resident_cumulative_strain_counts[0] + resident_cumulative_strain_counts[1]);

                delete[] health_state_change_time;
                delete[] current_health;
                delete[] current_location;
                delete[] current_strain;
                delete[] disease_profile_index_dict;
                delete[] resident_health_state_counts;
                delete[] resident_strain_counts;
                delete[] resident_cumulative_strain_counts;

            }
        }
    }

    delete[] resident;
    delete[] hospitals;
    delete[] cemeteries;
    delete[] pt_locations;
    delete[] pt_availability;
    delete[] initial_location;

    return 0;
}