from cell import Cell
from matplotlib import pyplot as plt
from matplotlib import use as mpl_use
from math import floor


def get_connection_factor(i, j):
    if 0 <= i <= 24 and 0 <= j <= 24:
        return [[0.6, 0.6, 0.6], [0.6, 0, 0.6], [0.6, 0.6, 0.6]]
    if 0 <= i <= 24 and 25 <= j <= 49:
        return [[1.0, 1.0, 1.0], [1.0, 0, 1.0], [1.0, 1.0, 1.0]]
    if 25 <= i <= 49 and 0 <= j <= 24:
        return [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    if 25 <= i <= 49 and 25 <= j <= 49:
        return [[0.3, 0.3, 0.3], [0.3, 0, 0.3], [0.3, 0.3, 0.3]]


class Space:
    def __init__(self, r: int, c: int, eps: float, virulence: float):
        # Defines an r x c grid of cells at time t=0
        self.r = r
        self.c = c
        self.t = 0
        self.population = 0

        self.susceptible = []
        self.infected = []
        self.recovered = []

        # Global parameters of the disease being modelled
        self.eps = eps
        self.virulence = virulence

        # Initialise 2D matrix of cells
        temp_m = [[0.5, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]]

        cells = []
        for i in range(r):
            row = []
            for j in range(c):
                # connection = get_connection_factor(i, j)
                connection = [[1.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]
                row.append(Cell([i, j], 100, connection, temp_m, 1.0, 0.0, 0.0))
                self.population += 100
            cells.append(row)
        cells[round(r / 2)-1][round(c / 2)-1].infected = [0.3]
        cells[round(r / 2)-1][round(c / 2)-1].susceptible = [0.7]

        self.cells = cells
        self.update_current_state()

    def __str__(self):
        return f"M: S:{round(self.susceptible[-1], 4)}, I:{round(self.infected[-1], 4)}, R:{round(self.recovered[-1], 4)}"

    # Returns the Moore neighbourhood of a given cell as a 2D array
    def get_moore_neighbourhood(self, coords: list[int]) -> list[list[any]]:
        neighbourhood = []

        for i in range(-1, 2):
            n = [None, None, None]
            row = coords[0] + i
            if not 0 <= row < self.r:
                neighbourhood.append(n)
                continue
            for j in range(-1, 2):
                column = coords[1] + j
                if not 0 <= column < self.c or (i == 0 and j == 0):
                    continue
                n[j + 1] = self.cells[row][column]
            neighbourhood.append(n)

        return neighbourhood

    # Returns the von Neumann neighbourhood of a given cell as a 2D array
    def get_vn_neighbourhood(self, coords: list[int]) -> list[list[any]]:
        neighbourhood = self.get_moore_neighbourhood(coords)

        neighbourhood[0][0] = None
        neighbourhood[0][2] = None
        neighbourhood[2][0] = None
        neighbourhood[2][2] = None
        return neighbourhood

    def neighbourhood_transition_term(self, neighbourhood: list[list[Cell]], cell: Cell) -> float:
        total = 0.0
        for row in range(3):
            for col in range(3):
                neighbour = neighbourhood[row][col]
                if neighbour is None:
                    continue

                c = cell.get_connection_factor(row, col)
                m = cell.get_movement_factor(row, col)
                total += (neighbour.population / cell.population) * c * m * self.virulence * neighbour.infected[self.t]

        return total

    def evolve(self):
        for r in self.cells:
            for cell in r:
                prev_i = cell.infected[self.t]
                prev_s = cell.susceptible[self.t]
                prev_r = cell.recovered[self.t]

                neighbourhood = self.get_moore_neighbourhood(cell.coords)
                n = self.neighbourhood_transition_term(neighbourhood, cell)

                s_to_i = self.virulence * prev_s * prev_i + prev_s * n
                if s_to_i > prev_s:
                    s_to_i = prev_s

                s = prev_s - s_to_i
                i = (1 - self.eps) * prev_i + s_to_i
                r = prev_r + self.eps * prev_i

                if not 0.99 <= s + i + r <= 1.01:
                    raise Exception("Something's gone wrong")

                cell.susceptible.append(s)
                cell.infected.append(i)
                cell.recovered.append(r)
                cell.discretise()

        self.update_current_state()
        self.t += 1

    def update_current_state(self):
        s = 0.0
        i = 0.0
        r = 0.0

        for row in range(self.r):
            for col in range(self.c):
                s += self.cells[row][col].susceptible[-1]
                i += self.cells[row][col].infected[-1]
                r += self.cells[row][col].recovered[-1]

        mean_s = s / (self.r * self.c)
        mean_i = i / (self.r * self.c)
        mean_r = 1 - mean_i - mean_s

        self.susceptible.append(mean_s)
        self.infected.append(mean_i)
        self.recovered.append(mean_r)

    def plot_sir_over_time(self):
        for i in range(self.t):
            print(
                f"T:{i}, S:{round(self.susceptible[i] * self.population)}, I:{round(self.infected[i] * self.population)}, "
                f"R:{round(self.recovered[i] * self.population)}")

        x = range(len(self.infected))

        mpl_use('MacOSX')
        plt.cla()
        plt.plot(x, self.infected, label="I")
        plt.plot(x, self.susceptible, label="S")
        plt.plot(x, self.recovered, label="R")
        plt.xlabel("t")
        plt.ylabel("Proportion of population")
        plt.legend()
        plt.show()

    def plot_state_at_times(self, times: list[int]):
        figure, axis = plt.subplots(2, 3)
        mpl_use('MacOSX')
        plt.cla()

        if len(times) > 6:
            return

        for t in times:
            i = []
            for r in range(self.r):
                row = []
                for c in range(self.c):
                    row.append(self.cells[r][c].discrete_infected[t])
                i.append(row)
            axis[floor(times.index(t) / 3), times.index(t) % 3].imshow(i)

        plt.show()
