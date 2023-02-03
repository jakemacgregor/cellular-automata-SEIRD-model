from space import Space

if __name__ == '__main__':
    columns = 50
    rows = 50
    eps = 0.4
    vir = 0.3
    iterations = 50

    if input("Do you want to specify parameters? (y/n)") == "y":
        columns = int(input("Number of columns (int):") or "50")
        rows = int(input("Number of rows (int):") or "50")
        eps = float(input("Epsilon value (float):") or "0.4")
        vir = float(input("Virulence (float):") or "0.3")
        iterations = int(input("Number of iterations (int):") or "50")

    space = Space(columns, rows, eps, vir)

    for i in range(iterations):
        space.evolve()

    space.plot_final_state()
