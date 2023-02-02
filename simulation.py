from space import Space

if __name__ == '__main__':
    space = Space(50, 50, 0.4, 0.3)

    for i in range(100):
        space.evolve()
        print(space.cells[24][24])

    space.print_plot_results()
