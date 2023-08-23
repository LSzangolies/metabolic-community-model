import mesa
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import pandas as pd

from mesa import Agent, Model
from mesa.time import BaseScheduler, RandomActivation
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector
from scipy.spatial.distance import euclidean

#function to calculate distance between cells in toroidal space:
def ToroidalDistance(p1, p2, model):
    x1, y1 = p1
    x2, y2 = p2

    dx = abs(x2 - x1)
    dy = abs(y2 - y1)

    if (dx > model.lengthx / 2):
        dx = model.lengthx - dx

    if (dy > model.lengthy / 2):
        dy = model.lengthy - dy

    return math.sqrt(dx * dx + dy * dy)

#individuals:
class mammal(Agent):
    def __init__(self, unique_id, x, y, spec, mass, real_mass, stor, age, model):
        super().__init__(unique_id, model)

        self.model = model
        self.mass = mass
        self.real_mass = real_mass
        self.storage = stor
        self.ss_use = False
        self.homequali = 0
        self.bold_prob = 1
        self.food_pref = "omni"
        self.for_type = "central"
        self.species = spec
        self.x = x
        self.y = y
        self.age = age
        self.mortbold = 0.333 / 10000  # background mortality
        self.move_factor = 3  # adaptation factor for movement distance (from data)
        self.shelter = 0  # individuals need no shelter,they use habitat edges
        self.hr_try_juv = 10  # attempts juveniles have to find a home range
        self.digestion = 0.5  # assimilation efficiency for herbivores
        self.sex = np.random.randint(0, 2)  # 0=male, 1=female
        self.growthcost = (7 + 6) * 1000 / 10  # costs for synthesizing flesh

        # allometrics
        self.update_allometrics(self.real_mass)
        self.maxhr_calc(self.real_mass)

        self.average_age = 1766.53 * (
                    self.mass ** 0.21)  # average lifespan: allometric formula of mammals after Hamilton et al 2011 [g, days]
        self.age_first_repro = 293.17 * (
                    self.mass ** 0.27)  # age of first reproduction: allometric fromula after Hamilton et al 2011
        self.gest_period = int(64.14 * (self.mass ** 0.24))  # gestation period allometric after Hamilton et al 2011
        self.lact_period = int(57.16 * (
                    self.mass ** 0.22))  # lactation period allometric after Hamilton et al 2011

        self.young = np.round(2.24 * self.mass ** -0.13)
        self.mass_neonate = 47.86 * self.mass ** 0.93
        self.mass_weaning = 295.12 * self.mass ** 0.91

        # initializations
        self.fmr = 0
        self.density = 0
        self.hrsize = 0
        self.living = True
        self.younggrowth = 0
        self.patchlist = list()
        self.feed_small = 0

        self.preg = 0
        self.pregcount = 0
        self.embryo_mass = 0
        self.real_young = 0

        if self.age > self.age_first_repro and self.sex == 1:
            self.preg = np.random.randint(0, 2)
            if self.preg == 1:
                self.real_young = self.young
                self.pregcount = np.random.randint(0, self.gest_period + self.lact_period)
                for i in range(min(self.pregcount, self.gest_period)):
                    self.embryo()
                if self.pregcount > self.gest_period:
                    for i in range(self.pregcount - self.gest_period):
                        self.juvenile()

    def feedrate_calc(self, m):
        self.feedrate = 25.6 * (m ** 0.737)  # daily feeding rate of mammals in dry g/day, mass in kg, after Nagy 2001

    def max_storage_calc(self, m):  # maximum energy storage capacity and foodshare
        self.max_storage = 3 * 294.8 * (m ** 1.19) + self.feedrate
        self.foodshare = (m / 0.001) ** -0.25

    def maxhr_calc(self, m):
        maxhr1 = 56.23 * (m ** 0.91)  # max hr for herbivores and omnivores, larger one is used
        maxhr2 = 47.863 * (m ** 1.18)
        self.maxhr = max(maxhr1, maxhr2) * 10000  # m2
        self.maxhr = math.sqrt(self.maxhr / math.pi)  # radius in m
        self.maxhr = self.maxhr / 10  # radius in patch length(10m)

    def update_allometrics(self, m):
        self.lococost = 10.7 * (m ** 0.68) + self.calc_postural(
            m)  # costs for mammals in J/m; mass in kg after Calder 1996
        self.lococost = self.lococost / 10000  # costs in g dry biomass/m; after Nagy'99 p.263 Buchmann 2011
        self.lococost = self.lococost * 10  # lococost in patchlength (= 10m)

        self.feedrate_calc(m)
        self.max_storage_calc(m)

    def calc_postural(self, m):  # postural locomotion costs
        posturaltime = (6.03 * (m ** 0.697) - 2.963 * (m ** 0.737))
        postural = posturaltime / (0.01398 * ((m * 1000) ** 0.217) * (np.e ** (-0.002 * ((np.log(m * 1000)) ** 2))))
        return (postural)

    def embryo(self):  # embryo growth curve
        A = np.exp(0.865 + 1.006 * np.log(self.mass_neonate))
        M0 = 0.0001
        b = np.log(A / M0)
        k = np.exp(0.627 - 0.905 * np.log(self.gest_period))
        self.embryo_mass = A * np.exp(-b * np.exp(-k * self.pregcount)) / 1000

    def juvenile(self):  # juvenile growth curve (during lactation)
        m = self.mass * 1000
        jm = self.embryo_mass * 1000
        factor = (m / jm) ** (1 / 3) - 1
        prod = factor * 3 / self.lact_period * np.log(
            (1 - ((self.mass_neonate / m) ** (1 / 3))) / (1 - ((self.mass_weaning / m) ** (1 / 3))))
        jm = prod * jm
        self.embryo_mass = self.embryo_mass + jm / 1000

    def growth(self):  # adult growth curve
        m = self.mass * 1000
        rm = self.real_mass * 1000
        factor = (m / rm) ** (1 / 3) - 1
        prod = factor * 3 / self.lact_period * np.log(
            (1 - ((self.mass_neonate / m) ** (1 / 3))) / (1 - ((self.mass_weaning / m) ** (1 / 3))))
        rm = prod * rm
        self.real_mass = min(self.mass, self.real_mass + rm / 1000)

        self.update_allometrics(self.real_mass)

    def ToroidalDistance(self, p1, p2):
        x1, y1 = p1
        x2, y2 = p2

        dx = abs(x2 - x1)
        dy = abs(y2 - y1)

        if (dx > self.model.lengthx / 2):
            dx = self.model.lengthx - dx

        if (dy > self.model.lengthy / 2):
            dy = self.model.lengthy - dy

        return math.sqrt(dx * dx + dy * dy)

    def find_hr(self):  # find suitable homerange for initial distribution
        success = 0
        tried = 0
        feed = 0
        while success == 0 and tried < 100:
            r_search = 0
            patch = [a[0][0] for a in self.model.grid.coord_iter() if a[0] != []]
            patch = np.random.choice([a for a in patch if isinstance(a, habitat) and a.home == True])
            tried = tried + 1
            self.x = patch.x
            self.y = patch.y
            self.model.grid.move_agent(self, (self.x, self.y))
            self.homequali = patch.patchquali
            feed = self.digestion * patch.patchfeed * self.foodshare
            if self.storage + feed >= self.max_storage:
                success = 1
                patch.spec_id = self.species
            self.patchlist = list([])
            r_search_local = 0
            while r_search_local < self.maxhr:
                local_patchlist = list()
                r_search_local = r_search_local + 1
                neighbors = patch.surroundings[r_search_local - 1]
                new_neighbors = list(set(neighbors) - set([item for subl in self.patchlist for item in subl]))
                self.patchlist.append(new_neighbors)
                if self.storage + feed < self.max_storage:
                    for patchy in new_neighbors:
                        if patchy.save == 1:
                            feed = feed + self.digestion * patchy.patchfeed * self.foodshare - 2 * self.lococost * self.ToroidalDistance(
                                (self.x, self.y), (patchy.x, patchy.y)) * self.move_factor
                        else:
                            feed = feed + (
                                        self.digestion * patchy.patchfeed * self.foodshare - 2 * self.lococost * self.ToroidalDistance(
                                    (self.x, self.y), (patchy.x, patchy.y)) * self.move_factor) * (1 - self.shelter)
                        r_search = r_search_local
            if self.storage + feed >= self.feedrate:
                success = 1
                self.storage = min(self.max_storage, feed - self.feedrate)
                patch.patchfeed = max(0, patch.patchfeed * (1 - self.foodshare))
                for neighbor in [item for subl in self.patchlist[0:(r_search - 1)] for item in subl]:
                    if neighbor.save == 1:
                        neighbor.patchfeed = max(0, neighbor.patchfeed * (1 - self.foodshare))
                    else:
                        neighbor.patchfeed = max(0, neighbor.patchfeed * (1 - self.foodshare * (1 - self.shelter)))
                    neighbor.spec_id = self.species

        if success == 1:
            self.hrsize = 2 * r_search
        if success == 0:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
            self.living = False

    def mother_feedrate(self):  # calculate energey need of reproducing females (dependent on juvenile mass and number)
        if self.age > 0:
            young_mass = 0
            if self.preg == 1:
                ym = self.embryo_mass
                if self.pregcount < self.gest_period:
                    self.embryo()
                    young_mass = self.real_young * self.embryo_mass
                else:
                    if self.pregcount < self.gest_period + self.lact_period:
                        self.juvenile()
                        young_mass = self.real_young * self.embryo_mass
                self.younggrowth = self.embryo_mass - ym
                self.feedrate_calc(young_mass + self.real_mass)
                self.max_storage_calc(young_mass + self.real_mass)

    def check_hr(self):  # daily foraging in home ranges
        self.mother_feedrate()
        if self.age > 0:
            self.lococost = self.lococost + (10.7 * (self.storage / 3930) ** 0.68 / 1000)
            success = 0
            r_search = 0
            this_cell = self.model.grid.get_cell_list_contents([(self.x, self.y)])
            patch = [obj for obj in this_cell if isinstance(obj, habitat)][0]
            feed = self.digestion * patch.patchfeed * self.foodshare
            self.fmr = self.fmr + 0.1 * (patch.patchfeed * self.foodshare)
            patch.patchfeed = max(0, patch.patchfeed * (1 - self.foodshare))
            if self.storage + feed >= self.max_storage:
                success = 1
                patch.spec_id = self.species
            while r_search < self.maxhr and self.storage + feed < self.max_storage:
                r_search = r_search + 1
                for neighbor in self.patchlist[(r_search - 1)]:
                    neighbor.spec_id = self.species
                    if neighbor.save == 1:
                        feed = feed + self.digestion * neighbor.patchfeed * self.foodshare - 2 * self.lococost * self.ToroidalDistance(
                            (self.x, self.y), (neighbor.x, neighbor.y)) * self.move_factor
                        self.fmr = self.fmr + 0.1 * (
                                    neighbor.patchfeed * self.foodshare) + 2 * self.lococost * self.ToroidalDistance(
                            (self.x, self.y), (neighbor.x, neighbor.y)) * self.move_factor
                        neighbor.patchfeed = neighbor.patchfeed * (1 - self.foodshare)
                    else:
                        feed = feed + (
                                    self.digestion * neighbor.patchfeed * self.foodshare - 2 * self.lococost * self.ToroidalDistance(
                                (self.x, self.y), (neighbor.x, neighbor.y)) * self.move_factor) * (1 - self.shelter)
                        self.fmr = self.fmr + (1 - self.shelter) * (0.1 * (
                                    neighbor.patchfeed * self.foodshare) + 2 * self.lococost * self.ToroidalDistance(
                            (self.x, self.y), (neighbor.x, neighbor.y)) * self.move_factor)
                        neighbor.patchfeed = neighbor.patchfeed * (1 - self.foodshare * (1 - self.shelter))
                    neighbor.patchfeed = max(0, neighbor.patchfeed)
            self.storage = min(self.max_storage, self.storage + feed)
            self.hrsize = 2 * r_search

    def maintenance(self):  # allocate ingested energy to all processes, die if not enough energy
        if self.age > 0:
            young_mass = 0
            synth = 0
            if self.preg == 1:
                young_mass = self.real_young * self.embryo_mass
                synth = self.real_young * self.younggrowth * self.growthcost
            old = self.feedrate
            self.feedrate_calc(young_mass + self.real_mass)
            act_feedrate = self.feedrate
            self.feedrate = old
            if self.storage > act_feedrate + synth:  # if there is enough food for maintenance
                m = self.real_mass
                self.growth()
                diff = self.real_mass - m
                diffcost = diff * self.growthcost
                self.storage = self.storage - act_feedrate - synth
                if self.storage > diffcost:  # if there is enough food to grow
                    self.storage = self.storage - diffcost
                    self.feedrate_calc(self.real_mass + young_mass)
                    act_feedrate = self.feedrate
                    self.feedrate_calc(self.real_mass)
                    self.fmr = self.fmr + act_feedrate + diffcost + synth
                else:
                    self.real_mass = m
                    self.update_allometrics(m)
                    self.feedrate_calc(self.real_mass + young_mass)
                    act_feedrate = self.feedrate
                    self.feedrate_calc(self.real_mass)
                    self.fmr = self.fmr + act_feedrate + synth
            else:
                if self.preg == 1:
                    if self.pregcount > self.gest_period:  # if there is not enough food and female is in lactation loose one juvenile after the other until enough
                        while (self.storage < self.feedrate + (
                                self.real_young * self.younggrowth * self.growthcost) and self.real_young > 0):
                            self.real_young = self.real_young - 1
                            young_mass = self.real_young * self.embryo_mass
                            self.feedrate_calc(young_mass + self.real_mass)
                        if self.storage > self.feedrate + (self.real_young * self.younggrowth * self.growthcost):
                            self.storage = self.storage - self.feedrate - (
                                        self.real_young * self.younggrowth * self.growthcost)
                            self.fmr = self.fmr + act_feedrate + (self.real_young * self.younggrowth * self.growthcost)
                            if self.real_young == 0:
                                self.preg = 0
                                self.pregcount = 0
                        else:  # if after loss of all juveniles there is still not enough: die
                            self.preg = 0
                            self.pregcount = 0
                            self.storage = 0
                            self.model.grid.remove_agent(self)
                            self.model.schedule.remove(self)
                            self.living = False
                    else:  # if not enough food and female is in gestation perios, loose pregnancy
                        self.preg = 0
                        self.pregcount = 0
                        self.feedrate_calc(self.real_mass)
                        if self.storage > self.feedrate:
                            self.storage = self.storage - self.feedrate
                            self.fmr = self.fmr + self.feedrate
                        else:  # if after loss of pregnancy there is still not enough: die
                            self.model.grid.remove_agent(self)
                            self.model.schedule.remove(self)
                            self.living = False
                else:  # if not enough food and not pregnant direct death
                    foodmort = (self.storage / self.feedrate)
                    self.storage = 0
                    self.model.grid.remove_agent(self)
                    self.model.schedule.remove(self)
                    self.living = False

    def offspring(self):  # initialization of new individuals after lactation
        stor = self.storage / (self.real_young + 1)
        real_mass = self.mass_weaning / 1000
        for i in range(int(self.real_young)):
            mass = (self.species + 1) * 1 / self.model.num_species * self.model.maxmass
            mass = np.random.normal(mass, 0.2 * mass)
            if mass > 0:
                off = mammal(self.model.next_id(), self.x, self.y, self.species, mass, real_mass, stor, 0, self.model)
                self.model.grid.place_agent(off, (off.x, off.y))
                off.maxhr_calc(off.mass)
                off.find_hr_offspring()
        self.storage = stor
        self.feedrate_calc(self.real_mass)
        self.max_storage_calc(self.real_mass)

    def get_food(self, neigh):  # function for food ingestion in juvenile home range search
        if neigh.save == 1:
            return self.digestion * neigh.patchfeed * self.foodshare - 2 * self.lococost * self.ToroidalDistance(
                (self.x, self.y), (neigh.x, neigh.y)) * self.move_factor
        else:
            return self.digestion * (
                        neigh.patchfeed * self.foodshare - 2 * self.lococost * self.ToroidalDistance((self.x, self.y), (
                neigh.x, neigh.y)) * self.move_factor) * (1 - self.shelter)

    def available_hr_offspring(self, patch):  # test wether there is a suitable home range for a juvenile at all
        r_search = math.ceil(self.maxhr)
        possible_feed = 0
        neighbors = patch.surroundings[int(self.maxhr)]
        (x, y) = (self.x, self.y)
        self.model.grid.move_agent(self, (patch.x, patch.y))
        possible_feed = map(self.get_food, neighbors)
        self.model.grid.move_agent(self, (x, y))
        possible_feed = sum(np.maximum(0, list(possible_feed)))
        return possible_feed > self.feedrate

    def reduce_food(self, neigh):  # function to reduce the resources after juvenile foraging
        if neigh.save == 1:
            neigh.patchfeed = neigh.patchfeed * (1 - self.foodshare)
        else:
            neigh.patchfeed = neigh.patchfeed * (1 - self.foodshare * (1 - self.shelter))
        neigh.spec_id = self.species

    def find_hr_offspring(self):  # actual search of juveniles for a home range
        self.foodshare = (self.mass / 0.001) ** (-0.25)  # allometric share of available food see Buchmann 2011
        max_nat_disp = 3.31 * (
                    self.mass ** 0.65)  # allometric maximum natal dispersal after Sutherland et al (2000) in km
        success = 0
        tried = 0
        feed = 0
        dispersal = self.model.grid.get_neighbors((self.x, self.y), moore=True, radius=int(100 * max_nat_disp))
        neighbors = [obj for obj in dispersal if
                     isinstance(obj, habitat) and self.ToroidalDistance((self.x, self.y), (obj.x, obj.y)) <= int(
                         100 * max_nat_disp) and obj.home == True]
        neighbors = np.random.choice(neighbors, min(len(neighbors), self.hr_try_juv))
        possible_hr = list(map(self.available_hr_offspring, neighbors))
        for patch in np.where(possible_hr)[0]:
            neighbor = neighbors[patch]
            r_search = 0
            self.x = neighbor.x
            self.y = neighbor.y
            self.model.grid.move_agent(self, (self.x, self.y))
            feed = self.digestion * neighbor.patchfeed * self.foodshare
            if self.storage + feed >= self.max_storage:
                success = 1
                neighbor.spec_id = self.species
            self.patchlist = list([])
            r_search_local = 0
            while r_search_local < self.maxhr:
                r_search_local = r_search_local + 1
                neighbors_here = neighbor.surroundings[r_search_local - 1]
                new_neighbors = list(set(neighbors_here) - set([item for subl in self.patchlist for item in subl]))
                if self.storage + feed < self.max_storage:
                    feed = feed + sum(map(self.get_food, new_neighbors))
                    r_search = r_search_local
                self.patchlist.append(new_neighbors)

            if self.storage + feed >= self.feedrate:
                success = 1
                self.storage = min(self.max_storage, self.storage + feed - self.feedrate)
                neighbor.patchfeed = neighbor.patchfeed * (1 - self.foodshare)
                list(map(self.reduce_food, [item for subl in self.patchlist[0:r_search] for item in subl]))
                break

        if success == 1:
            self.hrsize = 2 * r_search
            self.homequali = neighbor.patchquali
            self.model.schedule.add(self)
        if success == 0:
            self.model.grid.remove_agent(self)

    def mort(self):  # mortality due to age and background mortality
        if self.shelter == 0 and np.random.uniform(0, 1) < self.mortbold:
            self.model.grid.remove_agent(self)
            self.model.schedule.remove(self)
        else:
            if self.age > np.random.normal(self.average_age, 0.1 * self.average_age):
                self.model.grid.remove_agent(self)
                self.model.schedule.remove(self)

    def density_fact(self):  # calculate density of conspecifics around an individual for foraging order
        radius = self.model.grid.get_neighbors((self.x, self.y), moore=True, radius=int(self.maxhr))
        neighbors_own = [obj for obj in radius if
                         isinstance(obj, mammal) and self.ToroidalDistance((self.x, self.y), (obj.x, obj.y)) <= int(
                             self.maxhr) and obj.species == self.species]
        neighbors_all = [obj for obj in radius if
                         isinstance(obj, mammal) and self.ToroidalDistance((self.x, self.y), (obj.x, obj.y)) <= int(
                             self.maxhr)]
        if len(neighbors_all) > 0:
            self.density = len(neighbors_own) / len(neighbors_all)
        else:
            self.density = 0

    def step(self):  # daily timestep of an individual
        self.age = self.age + 1
        self.intake = 0
        self.fmr = 0
        self.density_fact()
        self.check_hr()
        self.maintenance()
        if self.living:
            if self.sex == 1:
                if self.preg == 1:
                    self.pregcount = self.pregcount + 1
                else:
                    if self.real_mass > 0.95 * self.mass:
                        self.preg = 1
                        self.pregcount = 0
                        self.real_young = self.young
                        self.embryo_mass = 0
                if self.pregcount > self.gest_period + self.lact_period:
                    self.offspring()
                    self.pregcount = 0
                    self.preg = 0
            self.mort()

#habitat cells:
class habitat(Agent):
    def __init__(self, unique_id, x, y, home, quali, model):
        super().__init__(unique_id, model)
        self.x = x
        self.y = y
        self.feed_struct = 16 * 1.7127
        self.spec_id = -1
        self.eaten = 0
        self.save = 0
        self.patchfeed = np.random.normal(self.feed_struct, 0.7)
        self.home = home
        self.patchquali = quali
        self.surroundings = [[]] * 20

    def step(self):  # daily timestep of the landscape - resources renew
        self.patchfeed = np.random.normal(self.feed_struct, 0.7)
        self.eaten = 0


class OrderActivation(BaseScheduler):  # scheduler defining the order of daily foraging of individuals
    def step(self):
        mass_order = 0.2
        maxmass = 0.1
        orders = np.array(
            [a.mass * (1 - a.density) for a in self.agents])  # there is a slight ordering after mass and density factor
        if self.get_agent_count() > 0:

            zuf = np.random.uniform(0, 1, self.get_agent_count())  # most individuals feed in random order
            orders[zuf > mass_order] = np.random.uniform(0, maxmass, np.sum(zuf > mass_order))

            for agent in np.array(self.agents)[np.argsort(-orders)]:
                agent.step()

class IBC(Model):  # actual full model definition
    def __init__(self, N, lengthx, lengthy, cover_p):
        self.individuals = N
        self.num_species = 10
        self.cover = cover_p * lengthx * lengthy  # habitat cover * # of patches (100 * 100)
        self.feed_matrix = 0.0  # no food prod in matrix
        self.mass_order = 0.2  # individuals that are ordered acording to body mass and density, remaining turtles are ordered randomly
        self.maxmass = 0.1
        self.lengthx = lengthx
        self.lengthy = lengthy
        self.current_id = 0
        self.grid = MultiGrid(lengthx, lengthy, True)
        self.schedule_patch = RandomActivation(self)
        self.schedule = OrderActivation(self)
        self.datacollector = DataCollector(model_reporters={
            "species_num": lambda m: len(np.unique([agent.species for agent in m.schedule.agent_buffer()])),
            "species 1": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 0]),
            "species 2": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 1]),
            "species 3": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 2]),
            "species 4": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 3]),
            "species 5": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 4]),
            "species 6": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 5]),
            "species 7": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 6]),
            "species 8": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 7]),
            "species 9": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 8]),
            "species 10": lambda m: len([agent for agent in m.schedule.agent_buffer() if agent.species == 9]),
            "species 1 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 0 and agent.age > 1]),
            "species 2 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 1 and agent.age > 1]),
            "species 3 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 2 and agent.age > 1]),
            "species 4 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 3 and agent.age > 1]),
            "species 5 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 4 and agent.age > 1]),
            "species 6 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 5 and agent.age > 1]),
            "species 7 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 6 and agent.age > 1]),
            "species 8 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 7 and agent.age > 1]),
            "species 9 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 8 and agent.age > 1]),
            "species 10 intake": lambda m: np.mean(
                [agent.intake for agent in m.schedule.agent_buffer() if agent.species == 9 and agent.age > 1]),
            "species 1 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 0 and agent.age > 1]),
            "species 2 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 1 and agent.age > 1]),
            "species 3 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 2 and agent.age > 1]),
            "species 4 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 3 and agent.age > 1]),
            "species 5 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 4 and agent.age > 1]),
            "species 6 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 5 and agent.age > 1]),
            "species 7 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 6 and agent.age > 1]),
            "species 8 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 7 and agent.age > 1]),
            "species 9 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 8 and agent.age > 1]),
            "species 10 hrsize": lambda m: np.mean(
                [agent.hrsize for agent in m.schedule.agent_buffer() if agent.species == 9 and agent.age > 1]),
            "species 1 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 0 and agent.age > 1]),
            "species 2 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 1 and agent.age > 1]),
            "species 3 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 2 and agent.age > 1]),
            "species 4 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 3 and agent.age > 1]),
            "species 5 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 4 and agent.age > 1]),
            "species 6 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 5 and agent.age > 1]),
            "species 7 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 6 and agent.age > 1]),
            "species 8 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 7 and agent.age > 1]),
            "species 9 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 8 and agent.age > 1]),
            "species 10 fmr": lambda m: np.mean(
                [agent.fmr for agent in m.schedule.agent_buffer() if agent.species == 9 and agent.age > 1])
            })  # output parameters

    def setup(self, fragmentation):
        self.setup_patches_algo(fragmentation)
        self.create_agents()

    def create_agents(self):  # initialization of individuals with their mass distribution
        zuf = np.random.uniform(0, 1, self.individuals)
        masses = (np.arange(1, self.num_species + 1) * 1 / self.num_species * self.maxmass) ** (-1.5)
        mass_sum = np.sum(masses)
        masses2 = np.zeros((10))
        for i in range(10):
            masses2[i] = np.sum(masses[0:(i + 1)]) / mass_sum
        num = 0
        for i in zuf:
            spec = np.sum(i > masses2)
            mass = (spec + 1) * 1 / self.num_species * self.maxmass
            mass = np.random.normal(mass, 0.2 * mass)
            if mass >= 0:
                age = 180 + np.random.randint(0, 180)
                x = self.random.randrange(self.grid.width)
                y = self.random.randrange(self.grid.height)
                a = mammal(num, x, y, spec, mass, mass, 0, age, self)
                self.schedule.add(a)
                self.grid.place_agent(a, (x, y))
                a.find_hr()
                num = num + 1
            else:
                self.dead = self.dead + 1
        self.current_id = self.individuals

    def setup_patches_algo(self,
                           clump):  # initialization of the landscape with habitat distributed according to the fragmentation level
        cov = 0
        random_order_agent = np.zeros((100 * 100))
        random_order_coord_x = np.zeros((100 * 100))
        random_order_coord_y = np.zeros((100 * 100))
        while cov < self.cover:
            positions = np.arange(100 * 100)
            np.random.shuffle(positions)
            i = 0
            for agent, x, y in self.grid.coord_iter():
                pos = positions[i]
                random_order_agent[pos] = len(agent)
                random_order_coord_x[pos] = int(x)
                random_order_coord_y[pos] = int(y)
                i = i + 1
            for i in range(100 * 100):
                agent = random_order_agent[i]
                x = int(random_order_coord_x[i])
                y = int(random_order_coord_y[i])
                if agent > 0:
                    neigh = self.grid.get_neighborhood((x, y), moore=True)
                    np.random.shuffle(neigh)
                    placed = False
                    for j in neigh:
                        if placed == False and (cov < self.cover) and self.grid.is_cell_empty(j):
                            patch = habitat(int(str(j[0]) + "0000" + str(j[1])), j[0], j[1], True, 1, self)
                            self.grid.place_agent(patch, (patch.x, patch.y))
                            self.schedule_patch.add(patch)
                            random_order_agent[
                                np.where(np.logical_and(random_order_coord_x == j[0], random_order_coord_y == j[1]))[0][
                                    0]] = 1
                            cov = cov + 1
                            placed = True
                else:
                    if np.random.uniform(0, 1) < (1 - clump) and (cov < self.cover):
                        patch = habitat(int(str(x) + "0000" + str(y)), x, y, True, 1, self)
                        self.grid.place_agent(patch, (patch.x, patch.y))
                        self.schedule_patch.add(patch)
                        cov = cov + 1

        # surroundings definition
        for agent in np.array(self.schedule_patch.agents):
            neighbors = self.grid.get_neighbors((agent.x, agent.y), moore=True)
            if np.sum([isinstance(obj, habitat) for obj in neighbors]) > 6:
                agent.save = 1
            for i in range(1, 21):
                agent.surroundings[i - 1] = list([obj for obj in
                                                  self.grid.get_neighbors((agent.x, agent.y), moore=True, radius=i,
                                                                          include_center=True) if
                                                  ToroidalDistance((agent.x, agent.y), (obj.x, obj.y), self) <= i])

    def step(self):  # daily step
        self.schedule_patch.step()
        self.schedule.step()
        self.datacollector.collect(self)

#simulating scenarios:
for frag in [0.9,0.99,0.999,0.9999]:
    for rep in range(10):
        model = IBC(1000, 100, 100, 0.05)
        model.setup(fragmentation=frag)
        for i in range(3650):
            model.step()
        result=model.datacollector.get_model_vars_dataframe()
        result.to_csv("Fragmentation"+str(frag)+"_Repetition"+str(rep)+".csv")
