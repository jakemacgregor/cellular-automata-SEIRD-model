from space import Space, plot_vaccination_results, write_to_csv
from matplotlib import pyplot as plt
from copy import deepcopy as copy

if __name__ == '__main__':
    columns = 50
    rows = 50
    sigma = 0.6
    eps = 0.4
    vir = 0.6
    iterations = 50
    output_timestamps = [0, 5, 10, 15, 20, 25]
    homogeneous_population = True
    constant_connection_factor = True
    constant_movement_factor = True
    start_in_center = True
    vaccination = False
    vaccination_time = 0
    vaccination_factors = [0.2, 0.3, 0.4]

    i_quarantine_factor = 0.0
    i_quarantine_trigger = 0.0
    e_quarantine_factor = 0.0
    e_quarantine_trigger = 0.0
    lockdown_trigger = 0.0
    unlock_trigger = 0.0

    # ENUM:
    uk_fast = False
    uk_slow = False

    if input("Use UK population data? (y/n)") == "y":
        columns = 109
        rows = 101
        uk_fast = True

        if input("Full 173 x 163 run? (y/n)") == "y":
            uk_slow = True
            uk_fast = False
            columns = 163
            rows = 173

    if input("Do you want to specify parameters? (y/n)") == "y":
        if not (uk_fast or uk_slow):
            columns = int(input("Number of columns (int):") or "50")
            rows = int(input("Number of rows (int):") or "50")
        sigma = float(input("Sigma value (float):") or "0.6")
        eps = float(input("Epsilon value (float):") or "0.4")
        vir = float(input("Virulence (float):") or "0.6")
        iterations = int(input("Number of iterations (int):") or "50")
        if not (uk_fast or uk_slow):
            start_in_center = not (input("Start infection in random location? (y/n)") == "y")
            homogeneous_population = not (input("Do you want inhomogeneous population distribution? (y/n)") == "y")
        constant_connection_factor = not (input("Do you want non-constant connections between cells? (y/n)") == "y")
        constant_movement_factor = not (input("Do you want non-constant movement between cells? (y/n)") == "y")

    if input("Do you want to simulate the effects of vaccination? (y/n)") == "y":
        vaccination = True
        vaccination_time = int(input("Timestep when vaccination begins:") or "16")

    if input("Do you want to simulate the effects of NPIs? (y/n)") == "y":
        i_quarantine_factor = float(input("Success rate of quarantining infected people (float):") or "0.98")
        i_quarantine_trigger = float(input("What % infected before quarantining infected people (float):") or "0.001")
        e_quarantine_factor = float(input("Success rate of asymptomatic quarantine (float):") or "0.2")
        e_quarantine_trigger = float(input("What % infected before asymptomatic quarantine (float):") or "0.003")
        lockdown_trigger = float(input("What % infected before restriction of movement (float):") or "0.05")
        unlock_trigger = float(input("What % infected before restriction of movement ends (float):") or "0.01")

    # Always create one space without vaccination
    spaces: list[Space] = [Space(rows, columns, sigma, eps, vir, 0, vaccination_time, i_quarantine_factor,
                                 i_quarantine_trigger, e_quarantine_factor, e_quarantine_trigger, lockdown_trigger,
                                 unlock_trigger, constant_connection_factor, homogeneous_population,
                                 constant_movement_factor, start_in_center, uk_fast, uk_slow)]

    for i in range(iterations):
        if i + 1 == vaccination_time and vaccination:
            a = [copy(spaces[0]) for i in range(3)]
            for j in range(3):
                a[j].set_vaccination_factor(vaccination_factors[j])
            spaces += a

        for space in spaces:
            space.evolve()

    write_to_csv(spaces[0])

    if input("Do you want to specify timestamps for cell space overview? (y/n)") == "y":
        output_timestamps = []
        print("Enter 6 integer timestamps:")
        for i in range(6):
            t = int(input())
            if not 0 <= t <= spaces[0].t:
                print("invalid timestamp")
                i -= 1
                continue
            output_timestamps.append(t)

    # Just plot the main graphs for the first space in the list to avoid getting too many figures to deal with
    spaces[0].plot_sir_over_time()
    spaces[0].plot_delta_sir_over_time()
    spaces[0].plot_infected_state_at_times(output_timestamps)
    spaces[0].plot_exposed_state_at_times(output_timestamps)
    if vaccination:
        plot_vaccination_results(spaces)
    plt.show()
