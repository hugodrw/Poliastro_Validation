from tqdm import tqdm
from ModSimulator import *
import numpy as np
from astropy import units as u


def make_data(action_list , initial_debris_list , starting_id , total_duration):
    """
    Input:
        action_list         : [(next_debris_norad_id , dt_given) , ...]
        initial_debris_list : [Debris]
        starting_id         : int
        total_duration      : int
    
    Output:
        dataframe : [frame_id | pos_otv | pos_debris_n , ...]
    """

    simulator = ModSimulator(starting_id , n_debris=len(initial_debris_list))

    for frame_id in tqdm(range(total_duration)):

        pass

    pass

def init_random_debris(n):
        """
        Output:
            list (norad_id , Orbit)
        """
        debris_list = []

        for norad_id in range(n):
            min_a = 6371 + 200
            max_a = 6371 + 2000
            a = np.random.uniform(min_a, max_a) * u.km
            ecc = 0 * u.one
            inc = np.random.uniform(0, 10) * u.deg
            raan = np.random.uniform(0, 10) * u.deg
            argp = 0 * u.deg
            nu = np.random.uniform(-180, 180) * u.deg

            debris = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)
            debris_list.append(Debris(poliastro_orbit=debris , norad_id=norad_id))

        return debris_list


if __name__ == "__main__":
    initial_debris_list = init_random_debris(n=10)
    simulator = ModSimulator(initial_debris_list , 1)
    a , b = simulator.strategy_1((9, 29))
    print(a)