# Cellular automata SEIRD model of infectious disease
*Created as my final year project*

This is a CA-based model capable of simulating the spread of an infectious disease through a population in an abstract space or in the UK, using UK Census data. To run:
1. Clone the repo
2. `pip install -r requirements.txt`
3. `python simulation.py` \ `python3 simulation.py`
4. Answer CLI prompts to tweak parameters
5. Graphical output can be seen in the `img` directory and will pop up at the end of the program running
6. A CSV file with the state at each timestep will also be created at the end, in the `csv_results` directory

Parameters:
- Sigma = proportion of exposed people who become infected in a given timestep
- Epsilon = proportion of infected people who recover in a given timestep
- Virulence = probability that contact between a sick and healthy person results in infection
- Xi = natural rate of birth and death within the population
- Zeta = proportion of infected people who die from the disease in a given timestep
- Remaining parameters are explained by the CLI prompts
