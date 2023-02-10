from cell import Cell
from matplotlib import pyplot as plt
from matplotlib import use as mpl_use
from math import floor
from random import seed, random


def get_connection_factor(i: int, j: int, const: bool) -> list[list[float]]:
    if const:
        return [[1.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]
    if 0 <= i <= 24 and 0 <= j <= 24:
        return [[0.6, 0.6, 0.6], [0.6, 0, 0.6], [0.6, 0.6, 0.6]]
    if 0 <= i <= 24 and 25 <= j:
        return [[1.0, 1.0, 1.0], [1.0, 0, 1.0], [1.0, 1.0, 1.0]]
    if 25 <= i and 0 <= j <= 24:
        return [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    if 25 <= i and 25 <= j:
        return [[0.3, 0.3, 0.3], [0.3, 0, 0.3], [0.3, 0.3, 0.3]]


def get_movement_factor(const: bool) -> list[list[float]]:
    if const:
        return [[0.5, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]]
    else:
        seed(123)
        return [[random(), random(), random()], [random(), 0, random()], [random(), random(), random()]]


def get_population(j: int, const: bool) -> int:
    if const:
        return 100
    else:
        return round(pow(1.17, j) * 10)


class Space:
    def __init__(self, r: int, c: int, eps: float, virulence: float, vaccination_factor: float, vaccination_time: int,
                 const_connection: bool, const_population: bool, const_movement: bool, start_center: bool):
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
        self.vaccination_factor = vaccination_factor
        self.vaccination_time = vaccination_time

        # Initialise 2D matrix of cells, setting the central cell to have 30% infected population
        cells = [[Cell([i, j], get_population(j, const_population), get_connection_factor(i, j, const_connection),
                       get_movement_factor(const_movement), 1.0, 0.0, 0.0) for j in range(c)] for i in range(r)]

        for row in cells:
            for cell in row:
                self.population += cell.population

        self.cells = cells
        self.start_infection(r, c, start_center)
        self.update_current_state()

    def __str__(self):
        return f"M: S:{round(self.susceptible[-1], 4)}, I:{round(self.infected[-1], 4)}, R:{round(self.recovered[-1], 4)}"

    def start_infection(self, r: int, c: int, center: bool) -> None:
        if center:
            i = round(r / 2) - 1
            j = round(c / 2) - 1
            self.cells[i][j].infected = [0.3]
            self.cells[i][j].susceptible = [0.7]
        else:
            seed(412)
            i = round(random() * r)
            j = round(random() * c)
            self.cells[i][j].infected = [0.3]
            self.cells[i][j].susceptible = [0.7]

    def set_vaccination_factor(self, factor: float) -> None:
        self.vaccination_factor = factor

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

    def evolve(self) -> None:
        for r in self.cells:
            for cell in r:
                prev_i = cell.infected[self.t]
                prev_s = cell.susceptible[self.t]
                prev_r = cell.recovered[self.t]

                neighbourhood = self.get_moore_neighbourhood(cell.coords)
                n = self.neighbourhood_transition_term(neighbourhood, cell)

                # Assume people are infected over being vaccinated as there may be some overlap
                s_to_i = self.virulence * prev_s * prev_i + prev_s * n
                if s_to_i > prev_s:
                    s_to_i = prev_s

                s_to_v = 0
                if self.t + 1 >= self.vaccination_time:
                    s_to_v = self.vaccination_factor * prev_s
                    if s_to_v > prev_s - s_to_i:
                        s_to_v = prev_s - s_to_i

                s = prev_s - s_to_v - s_to_i
                i = (1 - self.eps) * prev_i + s_to_i
                r = prev_r + self.eps * prev_i + s_to_v

                cell.susceptible.append(s)
                cell.infected.append(i)
                cell.recovered.append(r)
                cell.discretise()

        self.update_current_state()
        self.t += 1

    def update_current_state(self) -> None:
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

        self.susceptible.append(round(mean_s * self.population))
        self.infected.append(round(mean_i * self.population))
        self.recovered.append(round(mean_r * self.population))

    def plot_sir_over_time(self) -> None:
        x = range(len(self.infected))

        mpl_use('MacOSX')
        plt.figure()
        plt.plot(x, self.infected, label="I")
        plt.plot(x, self.susceptible, label="S")
        plt.plot(x, self.recovered, label="R")
        plt.xlabel("t")
        plt.ylabel("Number of people")
        plt.legend()
        plt.show()

    def plot_state_at_times(self, times: list[int]) -> None:
        figure, axis = plt.subplots(2, 3)
        mpl_use('MacOSX')

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


def plot_vaccination_results(s: list[Space]) -> None:
    plt.figure()
    x = range(len(s[0].infected))
    mpl_use('MacOSX')
    plt.plot(x, s[0].infected, label="0")
    plt.plot(x, s[1].infected, label="0.2")
    plt.plot(x, s[2].infected, label="0.3")
    plt.plot(x, s[3].infected, label="0.4")
    plt.xlabel("t")
    plt.ylabel("Number of people")
    plt.legend()

    for space in s:
        max_infected = max(space.infected)
        print(
            f"Vaccination: {space.vaccination_factor}, max I: {max_infected}, at time T: {space.infected.index(max_infected)}, total infected: {sum(space.infected)}")
