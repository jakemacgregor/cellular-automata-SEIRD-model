from cell import Cell


class Space:
    def __init__(self, r, c):
        # Defines an r x c grid of cells
        self.r = r
        self.c = c

        self.cells = [[] * c]
