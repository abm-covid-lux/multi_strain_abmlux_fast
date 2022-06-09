import os.path as osp
import logging
import logging.config
import math
from tqdm import tqdm
from ms_abmlux.sim_factory import SimulationFactory

log = logging.getLogger()

def build_integer_agents(sim):
    """Extracts data from a simulation factory and converts to numpy arrays of integers"""

    log.info("Extracing clock and agents...")

    # Build new sim and extract required objects
    agents = sim.world.agents
    rescale_factor = int(1/sim.world.scale_factor)

    # Agents
    number_of_agents = len(agents)
    resident = []
    for n in range(number_of_agents):
        if agents[n].region == "Luxembourg":
            resident.append(1)
        else:
            resident.append(0)

    return number_of_agents, resident, rescale_factor

def build_integer_locations(sim):
    """Extracts data from a simulation factory and converts to numpy arrays of integers"""

    log.info("Extracting locations...")

    agents = sim.world.agents
    number_of_agents = len(agents)
    locations = sim.world.locations
    transport_model = sim.transport_model
    scale_factor = sim.world.scale_factor

    # Locations
    number_of_locations = len(locations)
    location_int_dict = {}
    for m in range(number_of_locations):
        location_int_dict[locations[m]] = m
    hospitals = []
    cemeteries = []
    pt_locations = []
    for m in range(number_of_locations):
        if locations[m].typ == "Hospital":
            hospitals.append(location_int_dict[locations[m]])
        if locations[m].typ == "Cemetery":
            cemeteries.append(location_int_dict[locations[m]])
        if locations[m].typ == "Public Transport":
            pt_locations.append(location_int_dict[locations[m]])
    num_hospitals = len(hospitals)
    num_cemeteries = len(cemeteries)
    num_pt_locations = len(pt_locations)
    initial_location = []
    for n in tqdm(range(number_of_agents)):
        initial_location.append(location_int_dict[agents[n].current_location])
    wkendday_av = transport_model.units_available_weekend_day
    wkday_av = transport_model.units_available_week_day
    pt_availability = wkendday_av + (wkday_av * 5) + wkendday_av
    pt_availability = [max(math.ceil(a * scale_factor), 1) for a in pt_availability]

    return number_of_locations, num_hospitals, hospitals, num_cemeteries, cemeteries,
           num_pt_locations, pt_locations, pt_availability, initial_location

def build_integer_activities(sim):
    """Extracts data from a simulation factory and converts to arrays of integers"""

    log.info("Extracting activities...")

    agents = sim.world.agents
    number_of_agents = len(agents)
    locations = sim.world.locations
    activity_model = sim.activity_model

    # Activities
    number_of_locations = len(locations)
    location_int_dict = {}
    for m in range(number_of_locations):
        location_int_dict[locations[m]] = m
    number_of_activities = len(activity_model.activities)
    locations_for_activity = []
    num_locations_for_activity = []
    weekly_routine = []
    max_locations_for_activity = max([len(agents[n].locations_for_activity(a)) for a in
                                      range(number_of_activities) for n in range(number_of_agents)])
    for n in tqdm(range(number_of_agents)):
        locations_n = []
        num_locations_n = []
        for a in range(number_of_activities):
            locations_n_a = [location_int_dict[l] for l in agents[n].locations_for_activity(a)]
            num_locations_n.append(len(locations_n_a))
            locations_n_a += [-1]*(max_locations_for_activity - len(locations_n_a))
            locations_n.append(locations_n_a)
        locations_for_activity.append(locations_n)
        num_locations_for_activity.append(num_locations_n)
        weekly_routine.append(activity_model.weeks_by_agent[agents[n]].weekly_routine)

    return weekly_routine, num_locations_for_activity, locations_for_activity

def build_integer_disease(sim):
    """Extracts data from a simulation factory and converts to arrays of integers"""

    log.info("Extracting disease...")

    agents = sim.world.agents
    number_of_agents = len(agents)
    locations = sim.world.locations
    number_of_locations = len(locations)
    disease_model = sim.disease_model

    # Disease
    number_of_strains = len(disease_model.strains)
    transmission_probabilities = []
    for i in range(number_of_strains):
        asymptomatic = disease_model.strains[i].transmission_probability['ASYMPTOMATIC']
        symptomatic  = disease_model.strains[i].transmission_probability['INFECTED']
        transmission_probabilities.append([asymptomatic, symptomatic])
    transmission_multiplier = []
    for m in range(number_of_locations):
        if locations[m].typ in disease_model.no_transmission_locations:
            multiplier = 0
        else:
            if locations[m].typ in disease_model.reduced_transmission_locations:
                multiplier = disease_model.reduced_transmission_factor
            else:
                multiplier = 1
        transmission_multiplier.append(multiplier)
    health_state_int_dict = {'SUSCEPTIBLE': 0,
                             'EXPOSED': 1,
                             'ASYMPTOMATIC': 2,
                             'PREINFECTED': 3,
                             'INFECTED': 4,
                             'HOSPITALIZING': 5,
                             'VENTILATING': 6,
                             'RECOVERED': 7,
                             'DEAD': 8}
    disease_profile_dict = []
    disease_durations_dict = []
    max_disease_profile_length = max([max([len(disease_model.disease_profile_dict[agent][strain])
                                     for agent in agents]) for strain in disease_model.strains])
    def durations_map(dur):
        """Converts int and NoneType durations to only ints"""
        if isinstance(dur, float):
            return int(dur)
        else:
            return -1
    for n in tqdm(range(number_of_agents)):
        profile = []
        durations = []
        for i in range(number_of_strains):
            profile_i = [health_state_int_dict[h] for h in 
                            disease_model.disease_profile_dict[agents[n]][disease_model.strains[i]]]
            profile_i += [-1]*(max_disease_profile_length-len(profile_i))
            profile.append(profile_i)
            durations_i = [durations_map(dur) for dur in
                          disease_model.disease_durations_dict[agents[n]][disease_model.strains[i]]]
            durations_i += [-1]*(max_disease_profile_length-len(durations_i))
            durations.append(durations_i)
        disease_profile_dict.append(profile)
        disease_durations_dict.append(durations)
        health = health_state_int_dict['SUSCEPTIBLE']
    num_initial_cases_by_strain = []
    for i in range(number_of_strains):
        num = disease_model.num_initial_cases[disease_model.strains[i]]
        num_initial_cases_by_strain.append(num)

    return transmission_probabilities, transmission_multiplier, disease_profile_dict,
           disease_durations_dict, num_initial_cases_by_strain

def extract_sim(SIM_FACTORY_FILENAME):
    if osp.isfile(SIM_FACTORY_FILENAME):
        sim_factory = SimulationFactory.from_file(SIM_FACTORY_FILENAME)
        logging.config.dictConfig(sim_factory.config['logging'])
        log.info("Existing factory loaded from %s", SIM_FACTORY_FILENAME)
        sim = sim_factory.new_sim(None)
    else:
        log.warning("State file not found")
    return sim
