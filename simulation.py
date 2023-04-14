from space import Space, plot_vaccination_results
from compartment import Compartment
from matplotlib import pyplot as plt
from copy import deepcopy as copy
from math import floor

if __name__ == '__main__':
    # Variables describing the model at a high level and whether to use UK dataset
    columns = 50
    rows = 50
    iterations = 50
    output_timestamps = []
    uk_data = False

    # Parameters of the disease
    sigma = 0.6
    eps = 0.4
    vir = 0.6
    xi = 0.0005
    zeta = 0.005

    # Distribution and movement of populations
    homogeneous_population = True
    constant_connection_factor = True
    constant_movement_factor = True
    start_in_center = True

    # Vaccination parameters
    vaccination = False
    vaccination_time = 0
    vaccination_factors = [0.2, 0.3, 0.4]

    # Non-pharmaceutical intervention parameters
    i_quarantine_factor = 0.0
    i_quarantine_trigger = 0.0
    e_quarantine_factor = 0.0
    e_quarantine_trigger = 0.0
    lockdown_trigger = 1.0
    unlock_trigger = 0.0

    # Set up model for UK dataset
    if input("Use UK population data? (y/n)") == "y":
        columns = 163
        rows = 173
        uk_data = True
        sigma = 1 / 5
        eps = 1 / 7
        vir = 0.4
        iterations = 100

    # Allow user override for model setup and population movement/distribution
    if input("Do you want to specify parameters? (y/n)") == "y":
        if not uk_data:
            columns = int(input("Number of columns (int):") or "50")
            rows = int(input("Number of rows (int):") or "50")
        sigma = float(input("Sigma value (float):") or "0.6")
        eps = float(input("Epsilon value (float):") or "0.4")
        vir = float(input("Virulence (float):") or "0.6")
        xi = float(input("Xi value (float):") or "0.0005")
        zeta = float(input("Zeta value (float):") or "0.005")
        iterations = int(input("Number of iterations (int):") or "50")
        if not uk_data:
            start_in_center = not (input("Start infection in random location? (y/n)") == "y")
            homogeneous_population = not (input("Do you want inhomogeneous population distribution? (y/n)") == "y")
        constant_connection_factor = not (input("Do you want non-constant connections between cells? (y/n)") == "y")
        constant_movement_factor = not (input("Do you want non-constant movement between cells? (y/n)") == "y")

    # Allow user to enable vaccination and set parameter
    if input("Do you want to simulate the effects of vaccination? (y/n)") == "y":
        vaccination = True
        vaccination_time = int(input("Timestep when vaccination begins:") or "16")

    # Allow setup of NPIs with user-specified parameters
    if input("Do you want to simulate the effects of NPIs? (y/n)") == "y":
        i_quarantine_factor = float(input("Success rate of quarantining infected people (float):") or "0.81")
        i_quarantine_trigger = float(input("What % infected before quarantining infected people (float):") or "0.001")
        e_quarantine_factor = float(input("Success rate of asymptomatic quarantine (float):") or "0.17")
        e_quarantine_trigger = float(input("What % infected before asymptomatic quarantine (float):") or "0.003")
        lockdown_trigger = float(input("What % infected before restriction of movement (float):") or "0.05")
        unlock_trigger = float(input("What % infected before restriction of movement ends (float):") or "0.005")

    # Create an initial space with the parameters following user input
    # spaces is a list as further spaces can be added for comparison
    spaces: list[Space] = [Space(rows, columns, sigma, eps, vir, xi, zeta, 0, vaccination_time, i_quarantine_factor,
                                 i_quarantine_trigger, e_quarantine_factor, e_quarantine_trigger, lockdown_trigger,
                                 unlock_trigger, constant_connection_factor, homogeneous_population,
                                 constant_movement_factor, start_in_center, uk_data)]

    # Main loop of the program: evolve the space the required number of iterations
    # If vaccination is true then create copies of the space with the different parameters
    for i in range(iterations):
        if i + 1 == vaccination_time and vaccination:
            a = [copy(spaces[0]) for i in range(3)]
            for j in range(3):
                a[j].set_vaccination_factor(vaccination_factors[j])
            spaces += a

        for space in spaces:
            space.evolve()

    # Results are written to a CSV
    if input("Save to CSV? (y/n)") == "y":
        spaces[0].write_to_csv()

    # Allow users specify different timesteps for the snapshots of the space - grid is designed for 6 such figures
    # Otherwise divide intervals by 6 and use these equally spaced intervals
    if input("Do you want to specify timestamps for cell space overview? (y/n)") == "y":
        print("Enter 6 integer timestamps:")
        for i in range(6):
            t = int(input())
            if not 0 <= t <= spaces[0].t:
                print("invalid timestamp")
                i -= 1
                continue
            output_timestamps.append(t)
    else:
        interval = floor(iterations / 6)
        for i in range(1, 7):
            output_timestamps.append(interval*i)

    # Plot figures for the initial space by default - this avoids having too many figures
    # Save plots to img folder - overwritten after each run, so make copies if required
    plt.rcParams['figure.figsize'] = [6.8, 4.8]
    spaces[0].plot_population_over_time(False, list(Compartment))
    plt.savefig("./img/SEIRD.png", bbox_inches='tight')

    spaces[0].plot_population_over_time(False, [Compartment.EXPOSED, Compartment.INFECTED, Compartment.DECEASED])
    plt.savefig("./img/EID.png", bbox_inches='tight')

    plt.rcParams['figure.figsize'] = [6, 3.6]
    spaces[0].plot_state_at_times(output_timestamps, Compartment.INFECTED)
    plt.savefig("./img/infected_time.png", bbox_inches='tight')

    spaces[0].plot_state_at_times(output_timestamps, Compartment.EXPOSED)
    plt.savefig("./img/exposed_time.png", bbox_inches='tight')

    # Separately plot the results of vaccination for the different spaces if vaccination has taken place
    if vaccination:
        plt.rcParams['figure.figsize'] = [6, 4]
        plot_vaccination_results(spaces)
        plt.savefig("./img/vaccination.png", bbox_inches='tight')
    plt.show()
