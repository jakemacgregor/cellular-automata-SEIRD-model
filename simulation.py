import numpy as np
from space import Space
from cell import Cell

if __name__ == '__main__':
    space = Space(10, 10, 0.1, 0.3)
    space.cells[5][5].infected = [0.1]
    space.cells[5][5].susceptible = [0.9]

    for i in range(10):
        space.evolve()
        print(space.cells[5][5])

    print(space)
