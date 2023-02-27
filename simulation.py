from space import Space, plot_vaccination_results, write_to_csv
from matplotlib import pyplot as plt
from copy import deepcopy as copy

if __name__ == '__main__':
    columns = 50
    rows = 50
    sigma = 0.8
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

    if input("Do you want to specify parameters? (y/n)") == "y":
        columns = int(input("Number of columns (int):") or "50")
        rows = int(input("Number of rows (int):") or "50")
        sigma = float(input("Sigma value (float):") or "0.8")
        eps = float(input("Epsilon value (float):") or "0.4")
        vir = float(input("Virulence (float):") or "0.6")
        iterations = int(input("Number of iterations (int):") or "50")
        homogeneous_population = not (input("Do you want inhomogeneous population distribution? (y/n)") == "y")
        constant_connection_factor = not (input("Do you want non-constant connections between cells? (y/n)") == "y")
        constant_movement_factor = not (input("Do you want non-constant movement between cells? (y/n)") == "y")
        start_in_center = not (input("Start infection in random location? (y/n)") == "y")

    if input("Do you want to simulate the effects of vaccination? (y/n)") == "y":
        vaccination = True
        vaccination_time = int(input("Timestep when vaccination begins:") or "16")

    # Always create one space without vaccination
    spaces: list[Space] = [Space(rows, columns, sigma, eps, vir, 0, vaccination_time, constant_connection_factor,
                                 homogeneous_population, constant_movement_factor, start_in_center)]

    for i in range(iterations):
        if i + 1 == vaccination_time and vaccination:
            a = [copy(spaces[0]) for i in range(3)]
            for j in range(3):
                a[j].set_vaccination_factor(vaccination_factors[j])
            spaces += a

        for space in spaces:
            space.evolve()

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
    spaces[0].plot_state_at_times(output_timestamps)
    if vaccination:
        plot_vaccination_results(spaces)
    plt.show()
    write_to_csv(spaces[0])
