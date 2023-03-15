class Cell:
    def __init__(self, coords: list[int], population: int, connection: list[list[float]], movement: list[list[float]],
                 susceptible: float, exposed: float, infected: float, recovered: float, deceased: float):
        # Co-ords are an array [i, j] where the cell is in position [i, j] in cell space
        self.coords = tuple(coords)
        self.population = population
        self.empty = population == 0


        # Nested array representing the neighbourhood for a cell of the form:
        # [[NW,N,NE],[W,X,E],[SW,S,SE]] (compass directions)
        # Each compass direction is replaced with a value for the connection or movement factor for that neighbour
        self.connection = connection
        self.movement = movement

        # Proportions of population of the cell who are susceptible, infected, or recovered
        # Discretised so that there is a finite number of states: (0.00, 0.01, 0.02, ..., 0.99, 1.00)
        # Stored as a list --> susceptible[t] = susceptible population at time step t
        if not self.empty:
            self.susceptible = [susceptible]
            self.exposed = [exposed]
            self.infected = [infected]
            self.recovered = [recovered]
            self.deceased = [deceased]
        if self.empty:
            self.susceptible = []
            self.exposed = []
            self.infected = []
            self.recovered = []
            self.deceased = []

        self.discrete_susceptible = []
        self.discrete_exposed = []
        self.discrete_infected = []
        self.discrete_recovered = []
        self.discrete_deceased = []

        self.discretise()

    def __str__(self):
        if self.empty:
            return f"Empty cell"
        return f"S: {self.susceptible[-1]}, E:{self.exposed[-1]} I: {self.infected[-1]}, R: {self.recovered[-1]}, " \
               f"D: {self.deceased[-1]}"

    def get_connection_factor(self, a: int, b: int) -> float:
        return self.connection[a][b]

    def get_movement_factor(self, a: int, b: int) -> float:
        return self.movement[a][b]

    def discretise(self):
        if self.empty:
            return

        self.discrete_susceptible.append(round(self.susceptible[-1] * 100) / 100)
        self.discrete_exposed.append(round(self.exposed[-1] * 100) / 100)
        self.discrete_infected.append(round(self.infected[-1] * 100) / 100)
        self.discrete_deceased.append(round(self.deceased[-1]*100) / 100)
        self.discrete_recovered.append(1 - self.discrete_susceptible[-1] - self.discrete_exposed[-1]
                                       - self.discrete_infected[-1] - self.discrete_deceased[-1])
        return
