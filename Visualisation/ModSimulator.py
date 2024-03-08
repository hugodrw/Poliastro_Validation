from astropy import units as u

from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter3D

import CustomManeuvres as CustomManeuvres

import copy
import numpy as np
import pandas as pd


class Debris:
    def __init__(self , poliastro_orbit , norad_id):
        self.poliastro_orbit = poliastro_orbit
        self.norad_id = norad_id



class ModSimulator:
    def __init__(self , debris_list , starting_index=1):
        # Initialise the debris dictionary and assign the otv to an Orbit
        self.debris_list = debris_list#self.init_random_debris(n=n_debris) 
        self.otv_orbit = copy.copy(self.debris_list[starting_index].poliastro_orbit)
        

    def simulate_action(self , action):
        """
        Input:
            action : (next_debris_norad_id , dt_given)
        Output:
            DV_required , DT_required

        Can expand to testing multiple strategies
        Updates all objects in simulator
        When propagating, add at the end the (dt from action - dt required)
        """
        DV_required , DT_required = self.strategy_1(action)
        return DV_required , DT_required


    def strategy_1(self , action , n_frames_per_second = 30):
        """
        Strategy 1 defined in transfer strategies slides
        1. Inc
        2. Raan
        3. Hohmann
        """
        
        column_names = ["otv"] + [f'debris{i+1}' for i in range(len(self.debris_list))]
        df = pd.DataFrame(columns=column_names)



        # Set the target from the action
        target_debris = self.debris_list[action[0]].poliastro_orbit

        # ---- Inclination change
        inc_change = CustomManeuvres.simple_inc_change(self.otv_orbit, target_debris , do_propagation=False)

        # Get the transfer time of the hoh_phas
        transfer_time = inc_change.get_total_time()

        # Propagate all debris to the end of the transfer

        n_sec = transfer_time.to(u.s).value
        n_frames = n_sec * n_frames_per_second
        step = 1e2
        dt = 1/30 * u.s * step
        
        for frame_id in range(0 , int(n_frames)): #n_frames

            self.otv_orbit = self.otv_orbit.propagate(dt)
            otv_pos , _ = self.otv_orbit.rv()
            otv_pos = otv_pos.to(u.km).value / 6371
            
            new_data = [otv_pos]
            for i , debris in enumerate(self.debris_list):
                self.debris_list[i].poliastro_orbit = debris.poliastro_orbit.propagate(dt)

                deb_i_pos , _ = self.debris_list[i].poliastro_orbit.rv()
                deb_i_pos = deb_i_pos.to(u.km).value / 6371
                new_data.append(deb_i_pos)

            df.loc[frame_id] = new_data
            
            


        # Apply the maneuver to the otv
        self.otv_orbit = self.otv_orbit.apply_maneuver(inc_change)
        


        # ---- Raan change
        target_debris = self.debris_list[action[0]].poliastro_orbit
        raan_change = CustomManeuvres.simple_raan_change(self.otv_orbit, target_debris)

        # Get the transfer time of the hoh_phas
        transfer_time = raan_change.get_total_time()

        # Propagate all debris to the end of the transfer
        for i , debris in enumerate(self.debris_list):
            self.debris_list[i].poliastro_orbit = debris.poliastro_orbit.propagate(transfer_time)
        
        # Apply the maneuver to the otv
        self.otv_orbit = self.otv_orbit.apply_maneuver(raan_change)
        

        # ---- Hohmann
        target_debris = self.debris_list[action[0]].poliastro_orbit
        hoh_change = CustomManeuvres.hohmann_with_phasing(self.otv_orbit, target_debris)

        # Get the transfer time of the hoh_phas
        transfer_time = hoh_change.get_total_time()

        # Propagate all debris to the end of the transfer
        for i , debris in enumerate(self.debris_list):
            self.debris_list[i].poliastro_orbit = debris.poliastro_orbit.propagate(transfer_time)
        
        # Apply the maneuver to the otv
        self.otv_orbit = self.otv_orbit.apply_maneuver(hoh_change)
        

        # Total resources used
        total_dv = hoh_change.get_total_cost() + raan_change.get_total_cost() + inc_change.get_total_cost()
        min_time = hoh_change.get_total_time() + raan_change.get_total_time() + inc_change.get_total_time()


        # Propagate with the extra time after the action
        extra_time = action[1] * u.day - min_time
        # print(extra_time)
        # if extra_time.value > 0:
        #     self.otv_orbit = self.otv_orbit.propagate(extra_time)
        #     for i , debris in enumerate(self.debris_list):
        #         self.debris_list[i].poliastro_orbit = debris.poliastro_orbit.propagate(extra_time)

        #return total_dv , min_time
        return df



    def init_random_debris(self , n):
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
    


