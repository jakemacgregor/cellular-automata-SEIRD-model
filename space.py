from cell import Cell


def discretise(n):
    return round(n * 100) / 100


class Space:
    def __init__(self, r, c, eps, virulence):
        # Defines an r x c grid of cells at time t=0
        self.r = r
        self.c = c
        self.t = 0

        # Global parameters of the disease being modelled
        self.eps = eps
        self.virulence = virulence

        temp_m_c = [[0.0, 0.0, 0.0] * 3]
        # Initialise 2D matrix of cells
        cells = []
        for i in range(r):
            row = []
            for j in range(c):
                row.append(Cell([i, j], 100, temp_m_c, temp_m_c, 1.0, 0.0, 0.0))
            cells.append(row)

        self.cells = cells

    # Returns a 1D array of neighbours of a given cell
    def get_neighbourhood(self, coords):
        neighbourhood = []

        for i in range(-1, 2):
            n = [None, None, None]
            row = coords[0] + i
            if not 0 <= row < self.r:
                neighbourhood = neighbourhood.append(n)
                continue
            for j in range(-1, 2):
                column = coords[1] + j
                if not 0 <= column < self.c:
                    continue
                n[j + 1] = self.cells[row][column]
            neighbourhood = neighbourhood.append(n)

        return neighbourhood

    def neighbourhood_transition_term(self, neighbourhood, cell):
        total = 0
        for r in range(3):
            for c in range(3):
                neighbour = neighbourhood[r][c]
                if neighbour is None:
                    continue

                c = cell.get_connection_factor(r, c)
                m = cell.get_movement_factor(r, c)
                total += (neighbour.population / cell.population) * c * m * self.virulence * neighbour.infected[self.t]

        return total

    def evolve(self):
        for r in self.cells:
            for cell in r:
                prev_i = cell.infected[self.t]
                prev_s = cell.susceptible[self.t]
                prev_r = cell.recovered[self.t]

                neighbourhood = self.get_neighbourhood(cell.coords)
                n = self.neighbourhood_transition_term(neighbourhood, cell)
                i = discretise((1 - self.eps) * prev_i + self.virulence * prev_s * prev_i + prev_s * n)
                s = discretise(prev_s - self.virulence * prev_s * prev_i - prev_s * n)
                r = discretise(prev_r + self.eps * prev_i)

                if not 0.99 < s + i + r < 1.01:
                    raise Exception("Something's gone wrong")

                cell.susceptible.append(s)
                cell.infected.append(i)
                cell.recovered.append(r)

        self.t += 1
