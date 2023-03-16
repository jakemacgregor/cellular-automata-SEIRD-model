from enum import Enum


class Compartment(Enum):
    """
    Enum assigning a numerical value to each compartment of the model
    """
    SUSCEPTIBLE = 1
    EXPOSED = 2
    INFECTED = 3
    RECOVERED = 4
    DECEASED = 5
