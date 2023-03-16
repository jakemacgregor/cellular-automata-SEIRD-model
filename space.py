import csv
import math

import numpy as np
from datetime import datetime
from cell import Cell
from matplotlib import pyplot as plt
from matplotlib import use as mpl_use
from math import floor
from random import random
from compartment import Compartment

uk_start_locations = [(160, 139), (137, 75), (142, 112), (85, 87), (154, 70), (96, 16), (96, 25), (164, 70), (100, 115),
                      (110, 81), (46, 26), (51, 63), (87, 69), (70, 79), (120, 60)]


def get_connection_factor(i: int, j: int, const: bool) -> list[list[float]]:
    """
    Return the connection factors for a particular cell based on predetermined connection factors for different 'zones'
    in the cell space
    :param i: row of cell
    :param j: column of cell
    :param const: determines whether a constant factor should be used, i.e. independent of position of cell
    :return: a 3x3 matrix representing connection factors between a cell and its neighbours
    """
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
    """
    Return the movement factors for a particular cell
    :param const: determines whether a constant factor should be returned or the factor should be random
    :return: a 3x3 matrix representing movement factors between a cell and its neighbours
    """
    if const:
        return [[0.5, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]]
    else:
        return [[random(), random(), random()], [random(), 0, random()], [random(), random(), random()]]


def get_population(j: int, const: bool) -> int:
    """
    Returns the population of a cell based on its column
    :param j: column of cell
    :param const: determines whether every cell has the same population or eastern cells have a higher population
    :return: the population of a cell
    """
    if const:
        return 100
    else:
        return round(pow(1.17, j) * 10)


def get_pop_uk(i: int, j: int, uk: np.ndarray, fast: bool) -> int:
    """
    When using UK data, returns the population for each cell. Cells are made by combining 1x1km squares from the raw
    data. The way in which these are combined depends on whether a less precise 'fast' run is taking place.
    :param i: row of cell
    :param j: column of cell
    :param uk: data representing UK population data
    :param fast: boolean value determining whether 12x6km or 7x4km cells will be used
    :return: the population for the cell
    """
    population = 0
    if fast:
        rows = 12
        columns = 6
    else:
        rows = 7
        columns = 4

    for r in range(i * rows, (i + 1) * rows):
        if r > 1210:
            continue
        for c in range(j * columns, (j + 1) * columns):
            if c > 651:
                continue
            population += uk[r][c]

    return population


class Space:
    """
    A class representing the cell space made up of cells.

    r : int
        number of rows of cells
    c : int
        number of columns of cells
    t : int
        current timestep
    population : int
        total population of cell space

    susceptible : list[int]
        susceptible[t] is the number of susceptible people at time t
    exposed : list[int]
        exposed[t] is the number of exposed people at time t
    infected : list[int]
        infected[t] is the number of infected people at time t
    recovered : list[int]
        recovered[t] is the number of recovered people at time t

    delta_susceptible : list[int]
        delta_susceptible[t] is the difference in number of susceptible people between times t and t-1
    delta_exposed : list[int]
        delta_exposed[t] is the difference in number of exposed people between times t and t-1
    delta_infected : list[int]
        delta_infected[t] is the difference in number of infected people between times t and t-1
    delta_recovered : list[int]
        delta_recovered[t] is the difference in number of recovered people between times t and t-1

    sigma : float
        The rate at which infected people transition to infected
    eps : float
        The rate at which infected people recover
    virulence : float
        How likely a disease is to spread between people
    xi : float
        The natural rate of birth and death
    zeta : float
        The rate of death for infected people

    vaccination_factor : float
        The proportion of people who take up the vaccine in a given timestep
    vaccination_time : int
        The number of timesteps before a vaccine is available

    i_quarantine_factor : float
        How effectively infected people are quarantined
    i_quarantine_trigger : float
        The proportion of the population who become infected before infected people are quarantined
    i_quarantining_active : bool
        Whether quarantining for infected people is currently active
    e_quarantine_factor : float
        How effectively exposed (asymptomatic) people are quarantined
    e_quarantine_trigger : float
        The proportion of the population who become infected before exposed (asymptomatic) people are quarantined
    e_quarantining_active : bool
        Whether quarantining for infected people is currently active
    lockdown_trigger : float
        The proportion of the population who become infected before a lockdown is initiated
    unlock_trigger : float
        The proportion of the population who are infected before a lockdown is ended
    lockdown_active : bool
        Whether a lockdown is currently happening

    nonempty_cells : int
    cells : list[list[Cell]]
    """

    def __init__(self, r: int, c: int, sigma: float, eps: float, virulence: float, xi: float, zeta: float,
                 vaccination_factor: float, vaccination_time: int, i_quarantine_factor: float,
                 i_quarantine_trigger: float, e_quarantine_factor: float, e_quarantine_trigger: float,
                 lockdown_trigger: float, unlock_trigger: float, const_connection: bool, const_population: bool,
                 const_movement: bool, start_center: bool, uk_fast: bool, uk_slow: bool):

        self.r = r
        self.c = c
        self.t = 0
        self.population = 0

        self.susceptible = []
        self.exposed = []
        self.infected = []
        self.recovered = []
        self.deceased = []

        self.delta_susceptible = [0]
        self.delta_exposed = [0]
        self.delta_infected = [0]
        self.delta_recovered = [0]
        self.delta_deceased = [0]

        self.sigma = sigma
        self.eps = eps
        self.virulence = virulence
        self.xi = xi
        self.zeta = zeta

        self.vaccination_factor = vaccination_factor
        self.vaccination_time = vaccination_time

        self.i_quarantine_factor = i_quarantine_factor
        self.i_quarantine_trigger = i_quarantine_trigger
        self.i_quarantining_active = False
        self.e_quarantine_factor = e_quarantine_factor
        self.e_quarantine_trigger = e_quarantine_trigger
        self.e_quarantining_active = False
        self.lockdown_trigger = lockdown_trigger
        self.unlock_trigger = unlock_trigger
        self.lockdown_active = False

        if uk_fast or uk_slow:
            data = np.loadtxt("UK_population.asc", skiprows=6)
            cells = [[Cell([i, j], get_pop_uk(i, j, data, uk_fast), get_connection_factor(i, j, const_connection),
                           get_movement_factor(const_movement), susceptible=1.0, exposed=0.0, infected=0.0,
                           recovered=0.0, deceased=0.0)
                      for j in range(c)] for i in range(r)]
        else:
            cells = [[Cell([i, j], get_population(j, const_population), get_connection_factor(i, j, const_connection),
                           get_movement_factor(const_movement), susceptible=1.0, exposed=0.0, infected=0.0,
                           recovered=0.0, deceased=0.0)
                      for j in range(c)] for i in range(r)]

        self.nonempty_cells = 0
        for row in cells:
            for cell in row:
                self.population += cell.population
                if not cell.empty:
                    self.nonempty_cells += 1

        self.cells = cells

        if uk_fast or uk_slow:
            self.start_infection_particular(uk_start_locations)
            # for i in range(15):
            #     self.start_infection_uk(r, c)
        else:
            self.start_infection(r, c, start_center)
        self.update_current_state()

    def __str__(self):
        return f"T:{self.t}, S:{round(self.susceptible[-1], 4)}, E:{round(self.exposed[-1], 4)}," \
               f" I:{round(self.infected[-1], 4)}, R:{round(self.recovered[-1], 4)}, D:{round(self.deceased[-1], 4)}"

    def start_infection(self, r: int, c: int, center: bool) -> None:
        """
        Begins the infection either in the center or a random position. Recursively calls itself until it finds a
        non-empty cell in which to start the infection. Sets 30% of the population to exposed.
        :param r: row of cell
        :param c: column of cell
        :param center: whether to start the infection in the center
        """
        if center:
            i = round(r / 2) - 1
            j = round(c / 2) - 1
        else:
            i = math.floor(random() * r)
            j = math.floor(random() * c)

        cell = self.cells[i][j]
        if cell.empty and not center:
            self.start_infection(r, c, center)
            return

        cell.susceptible = [0.7]
        cell.exposed = [0.3]

    def start_infection_particular(self, locations: list[tuple[int, int]]):
        for location in locations:
            cell = self.cells[location[0]][location[1]]
            cell.susceptible = [0.7]
            cell.exposed = [0.3]

    def start_infection_uk(self, r: int, c: int) -> None:
        """
        Similar to start_infection, but specifically for UK data. Recursively calls itself until it finds a cell with
        a population of at least 100. Sets 30% of the population to exposed.
        :param r: row of cell
        :param c: column of cell
        """
        i = math.floor(random() * r)
        j = math.floor(random() * c)

        cell = self.cells[i][j]
        if cell.empty or cell.population < 100 or not cell.susceptible[0] == 1.0:
            self.start_infection_uk(r, c)
            return

        cell.susceptible = [0.7]
        cell.exposed = [0.3]

    def set_vaccination_factor(self, factor: float) -> None:
        self.vaccination_factor = factor

    # Returns the Moore neighbourhood of a given cell as a 2D array
    def get_moore_neighbourhood(self, coords: list[int]) -> list[list[any]]:
        """
        Returns the Moore neighbourhood of a given cell as a 2D array
        """
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
        """
        Returns the von Neumann neighbourhood of a given cell as a 2D array
        """
        neighbourhood = self.get_moore_neighbourhood(coords)

        neighbourhood[0][0] = None
        neighbourhood[0][2] = None
        neighbourhood[2][0] = None
        neighbourhood[2][2] = None
        return neighbourhood

    def neighbourhood_transition_term(self, neighbourhood: list[list[Cell]], cell: Cell) -> float:
        """
        Calculates the effect of a cell's neighbourhood on the number of people in a cell who are exposed to the
        disease. For full context see transition functions in documentation.
        """
        total = 0.0
        for row in range(3):
            for col in range(3):
                neighbour = neighbourhood[row][col]
                if neighbour is None:
                    continue

                if neighbour.empty:
                    continue

                if self.lockdown_active:
                    c = 0.1
                    m = 0.1
                else:
                    c = cell.get_connection_factor(row, col)
                    m = cell.get_movement_factor(row, col)

                infected = neighbour.infected[self.t]
                if self.i_quarantining_active:
                    infected = (1 - self.i_quarantine_factor) * infected

                exposed = neighbour.exposed[self.t]
                if self.e_quarantining_active:
                    exposed = (1 - self.e_quarantine_factor) * exposed

                total += (neighbour.population / cell.population) * c * m * self.virulence * \
                         (0.5 * exposed + infected)

        return total

    def evolve(self) -> None:
        """
        Applies the transition function to each cell in the cell space and stores the updated states of each. See the
        full transition functions in the documentation for a thorough explanation.
        """
        for r in self.cells:
            for cell in r:
                if cell.empty:
                    continue

                prev_s = cell.susceptible[self.t] + self.xi * (1 - cell.susceptible[self.t])
                prev_e = cell.exposed[self.t] - self.xi * cell.exposed[self.t]
                prev_i = cell.infected[self.t] - self.xi * cell.infected[self.t]
                prev_r = cell.recovered[self.t] - self.xi * cell.recovered[self.t]
                prev_d = cell.deceased[self.t]

                neighbourhood = self.get_moore_neighbourhood(cell.coords)
                n = self.neighbourhood_transition_term(neighbourhood, cell)

                infected_minus_quarantine = prev_i
                if self.i_quarantining_active:
                    infected_minus_quarantine = prev_i * (1 - self.i_quarantine_factor)

                exposed_minus_quarantine = prev_e
                if self.e_quarantining_active:
                    exposed_minus_quarantine = prev_e * (1 - self.e_quarantine_factor)

                # Assume people are infected over being vaccinated as there may be some overlap
                s_to_e = self.virulence * prev_s * (infected_minus_quarantine + 1 / 2 * exposed_minus_quarantine) \
                         + prev_s * n
                if s_to_e > prev_s:
                    s_to_e = prev_s

                s_to_v = 0
                if self.t + 1 >= self.vaccination_time:
                    s_to_v = self.vaccination_factor * prev_s
                    if s_to_v > prev_s - s_to_e:
                        s_to_v = prev_s - s_to_e

                s = prev_s - s_to_v - s_to_e
                e = (1 - self.sigma) * prev_e + s_to_e
                i = (1 - (self.eps + self.zeta)) * prev_i + self.sigma * prev_e
                r = prev_r + self.eps * prev_i + s_to_v
                d = prev_d + self.zeta * prev_i

                cell.susceptible.append(s)
                cell.exposed.append(e)
                cell.infected.append(i)
                cell.recovered.append(r)
                cell.deceased.append(d)
                cell.discretise()

        self.update_current_state()
        self.update_delta_values()
        self.t += 1

    def update_current_state(self) -> None:
        """
        Called at the end of an evolution, calculates the mean state of the cells and the approximate S, E, I, R, D
        populations of the cell space overall. Also determines if thresholds have been met to trigger NPIs
        """
        s = 0.0
        e = 0.0
        i = 0.0
        r = 0.0
        d = 0.0

        for row in range(self.r):
            for col in range(self.c):
                cell = self.cells[row][col]
                if cell.empty:
                    continue
                s += cell.susceptible[-1]
                e += cell.exposed[-1]
                i += cell.infected[-1]
                r += cell.recovered[-1]
                d += cell.deceased[-1]

        mean_s = s / self.nonempty_cells
        mean_e = e / self.nonempty_cells
        mean_i = i / self.nonempty_cells
        mean_d = d / self.nonempty_cells
        mean_r = 1 - mean_e - mean_i - mean_s - mean_d

        if mean_i >= self.i_quarantine_trigger:
            self.i_quarantining_active = True

        if mean_i >= self.e_quarantine_trigger:
            self.e_quarantining_active = True

        if mean_i >= self.lockdown_trigger:
            self.lockdown_active = True
        if mean_i <= self.unlock_trigger:
            self.lockdown_active = False

        self.susceptible.append(round(mean_s * self.population))
        self.exposed.append(round(mean_e * self.population))
        self.infected.append(round(mean_i * self.population))
        self.recovered.append(round(mean_r * self.population))
        self.deceased.append(round(mean_d * self.population))

    def update_delta_values(self) -> None:
        """
        Updates the delta values
        """
        s, e, i, r, d = self.susceptible[-1], self.exposed[-1], self.infected[-1], self.recovered[-1], self.deceased[-1]
        prev_s, prev_e, prev_i, prev_r, prev_d = self.susceptible[-2], self.exposed[-2], self.infected[-2], \
            self.recovered[-2], self.deceased[-2]

        self.delta_susceptible.append((s - prev_s))
        self.delta_exposed.append((e - prev_e))
        self.delta_infected.append((i - prev_i))
        self.delta_recovered.append((r - prev_r))
        self.delta_deceased.append((d - prev_d))

    def plot_population_over_time(self, delta: bool, compartments: list[Compartment]) -> None:
        """
        Plot chosen compartments over each time step.
        """
        x = range(len(self.infected))
        mpl_use('MacOSX')

        plt.figure()
        for compartment in compartments:
            match compartment:
                case Compartment.SUSCEPTIBLE:
                    data = self.susceptible
                    if delta:
                        data = self.delta_susceptible
                    plt.plot(x, data, c='limegreen', label="S")
                case Compartment.EXPOSED:
                    data = self.exposed
                    if delta:
                        data = self.delta_exposed
                    plt.plot(x, data, c='gold', label="E")
                case Compartment.INFECTED:
                    data = self.infected
                    if delta:
                        data = self.delta_infected
                    plt.plot(x, data, c='orangered', label="I")
                case Compartment.RECOVERED:
                    data = self.recovered
                    if delta:
                        data = self.delta_recovered
                    plt.plot(x, data, c='darkviolet', label="R")
                case Compartment.DECEASED:
                    data = self.deceased
                    if delta:
                        data = self.delta_deceased
                    plt.plot(x, data, c='dimgray', label="D")

        plt.xlabel("t")
        plt.ylabel("Number of people")
        plt.legend()
        plt.show()

    def plot_state_at_times(self, times: list[int], compartment: Compartment) -> None:
        """
        Plot a heatmap of the cell space at particular times for a particular compartment
        :param times: times at which to plot the heatmaps
        :param compartment: compartment to plot
        """
        mpl_use('MacOSX')

        figure, axis = plt.subplots(2, 3)
        if len(times) > 6:
            return

        max_value = 0.0
        for t in times:
            i = []
            for r in range(self.r):
                row = []
                for c in range(self.c):
                    cell = self.cells[r][c]
                    if cell.empty:
                        row.append(0)
                    else:
                        d = 0
                        match compartment:
                            case Compartment.SUSCEPTIBLE:
                                d = cell.discrete_susceptible[t]
                            case Compartment.EXPOSED:
                                d = cell.discrete_exposed[t]
                            case Compartment.INFECTED:
                                d = cell.discrete_infected[t]
                            case Compartment.RECOVERED:
                                d = cell.discrete_recovered[t]
                            case Compartment.DECEASED:
                                d = cell.discrete_deceased[t]
                        if d > max_value:
                            max_value = d
                        row.append(d)
                i.append(row)

            im = axis[floor(times.index(t) / 3), times.index(t) % 3].imshow(i, vmin=0, vmax=max_value, cmap='plasma')
            axis[floor(times.index(t) / 3), times.index(t) % 3].set_title(f"t = {t}")
            axis[floor(times.index(t) / 3), times.index(t) % 3].axis('off')

        figure.subplots_adjust(right=0.8)
        cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
        figure.colorbar(im, cax=cbar_ax)

    def write_to_csv(self) -> None:
        """
        Records the state of the cell space at each time step in a CSV file to make more detailed analysis easy
        """
        header = ['t', 's', 'ds', 'e', 'de', 'i', 'di', 'r', 'dr', 'd', 'dd']

        data = []
        for i in range(self.t):
            row = [f"{i}", f"{self.susceptible[i]}", f"{self.delta_susceptible[i]}", f"{self.exposed[i]}",
                   f"{self.delta_exposed[i]}", f"{self.infected[i]}", f"{self.delta_infected[i]}",
                   f"{self.recovered[i]}", f"{self.delta_recovered[i]}", f"{self.deceased[i]}",
                   f"{self.delta_deceased[i]}"]
            data.append(row)

        dt = datetime.now().strftime("%Y-%m-%d %H;%M;%S")
        with open(f"./csv_results/{dt}.csv", 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(data)


def plot_vaccination_results(s: list[Space]) -> None:
    """
    Shows how trends in infected populations differ with the different vaccination parameters
    """
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
            f"Vaccination: {space.vaccination_factor}, max I: {max_infected}, at time T: "
            f"{space.infected.index(max_infected)}, total infected: {sum(space.infected)}")
