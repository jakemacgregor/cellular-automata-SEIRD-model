from cell import Cell

neighbours = [["NW", "N", "NE"], ["W", "X", "E"], ["SW", "S", "SE"]]


class Space:
    def __init__(self, r, c):
        # Defines an r x c grid of cells
        self.r = r
        self.c = c

        cells = []
        for j in range(c):
            cells.append([])

        self.cells = cells

    def get_neighbourhood(self, coords):
        neighbourhood = {}

        for i in range(-1, 2):
            if not 0 <= coords[0] + i < self.r:
                continue
            for j in range(-1, 2):
                if not 0 <= coords[1] + j < self.c:
                    continue
                n = neighbours[i][j]
                neighbourhood[n] = [i, j]

        return neighbourhood
