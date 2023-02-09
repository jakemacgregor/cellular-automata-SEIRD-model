from space import Space
from matplotlib import pyplot as plt
from matplotlib import use as mpl_use

if __name__ == '__main__':
    columns = 50
    rows = 50
    eps = 0.4
    vir = 0.6
    iterations = 50
    output_timestamps = [0, 5, 10, 15, 20, 25]
    homogeneous_population = True
    constant_connection_factor = True
    constant_movement_factor = True
    start_in_center = True

    if input("Do you want to specify parameters? (y/n)") == "y":
        columns = int(input("Number of columns (int):") or "50")
        rows = int(input("Number of rows (int):") or "50")
        eps = float(input("Epsilon value (float):") or "0.4")
        vir = float(input("Virulence (float):") or "0.6")
        iterations = int(input("Number of iterations (int):") or "50")
        homogeneous_population = not (input("Do you want inhomogeneous population distribution? (y/n)") == "y")
        constant_connection_factor = not (input("Do you want non-constant connections between cells? (y/n)") == "y")
        constant_movement_factor = not (input("Do you want non-constant movement between cells? (y/n)") == "y")
        start_in_center = not (input("Start infection in random location? (y/n)") == "y")

    space = Space(columns, rows, eps, vir, 0, 16, homogeneous_population, constant_connection_factor,
                  constant_movement_factor, start_in_center)

    space_2 = Space(columns, rows, eps, vir, 0.2, 16, homogeneous_population, constant_connection_factor,
                    constant_movement_factor, start_in_center)

    space_3 = Space(columns, rows, eps, vir, 0.3, 16, homogeneous_population, constant_connection_factor,
                    constant_movement_factor, start_in_center)

    space_4 = Space(columns, rows, eps, vir, 0.4, 16, homogeneous_population, constant_connection_factor,
                    constant_movement_factor, start_in_center)

    for i in range(iterations):
        space.evolve()
        space_2.evolve()
        space_3.evolve()
        space_4.evolve()

    if input("Do you want to specify timestamps for graphical output? (y/n)") == "y":
        output_timestamps = []
        print("Enter 6 integer timestamps:")
        for i in range(6):
            t = int(input())
            if not 0 <= t <= space.t:
                print("invalid timestamp")
                i -= 1
                continue
            output_timestamps.append(t)

    # space.plot_sir_over_time()
    # space.plot_state_at_times(output_timestamps)
    x = range(len(space.infected))
    mpl_use('MacOSX')
    plt.cla()
    plt.plot(x, space.infected, label="0")
    plt.plot(x, space_2.infected, label="0.2")
    plt.plot(x, space_3.infected, label="0.4")
    plt.plot(x, space_4.infected, label="0.6")
    plt.xlabel("t")
    plt.ylabel("Number of people")
    plt.legend()
    plt.show()