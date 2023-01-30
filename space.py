from cell import Cell


class Space:
    def __init__(self, r, c):
        # Defines an r x c grid of cells
        self.r = r
        self.c = c

        self.cells = [[] * c]

    def get_neighbourhood(self, cell_coords):
        cell = self.cells[cell_coords[0][cell_coords[1]]]
        neighbours = [[], [], []]

        for i in range(-1, 2):
            if cell_coords[0] + i < 0 or cell_coords[0] + i >= self.r:
                continue
            for j in range(-1, 2):
                if cell_coords[0] + j < 0 or cell_coords[0] + j >= self.c:
                    continue
                neighbours[(i+1)[j+1]] = [cell_coords[0] + i, cell_coords[1] + j]
