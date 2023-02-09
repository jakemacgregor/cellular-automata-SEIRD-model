from space import Space

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

    space = Space(columns, rows, eps, vir, homogeneous_population, constant_connection_factor, constant_movement_factor,
                  start_in_center)

    for i in range(iterations):
        space.evolve()

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

    space.plot_sir_over_time()
    space.plot_state_at_times(output_timestamps)
