import numpy as np
from space import Space
from cell import Cell

if __name__ == '__main__':
    space = Space(10, 10, 0.1, 0.3)

    for i in range(100):
        space.evolve()
        print(space.cells[5][5])
        print(space)

    space.print_results()