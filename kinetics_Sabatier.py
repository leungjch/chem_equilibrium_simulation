from __future__ import division
import pygame, sys, pygame.mixer
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import pylab
import numpy as np
import time
import random
from pygame.locals import *
import pylab
plt.style.use("ggplot")
# Add 75 CO2 and 15 H2
WIDTH, HEIGHT = 1920//2,900
SPEED =  6
conc_CO2 = input("Enter the initial concentration of CO2\n"+">")

conc_H2 = input("Enter the initial concentration of H2\n"+">")

conc_CH4 = input("Enter the initial concentration of CH4\n"+">")

conc_H2O = input("Enter the initial concentration of H2O\n"+">")

add = 20
showGraph = False
# conc_CO2 =20
# conc_H2 = 20
# conc_CH4 =0
# conc_H2O =1
tickrate = 10
V = 1

n_CO2 = int(float(conc_CO2) * V)
n_H2 = int(float(conc_H2) * V)
n_CH4 = int(float(conc_CH4) * V)
n_H2O = int(float(conc_H2O)* V)

radius_CO2 = 330/265 * 10
radius_H2 = 289/265 * 10
radius_CH4 = 380/265 * 10
radius_H2O = 265/265 * 10

mass_CO2 = 44
mass_H2 = 2
mass_CH4 = 16
mass_H2O = 18


# Act energy
Ea_H2_H2 = 0
Ea_4H2_CO2 = 5
Ea_H2O_CO2 = 10
Ea_H2O_CO2 = 10

Collision_H2_CO2 = ["H2", "CO2"]
Collision_H2O_CH4 = ["CO2", "CH4"]


pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))

clock = pygame.time.Clock()



molecules_CO2 = []
molecules_H2 = []
molecules_CH4 = []
molecules_H2O = []
molecules_All = []

colour_CO2 = (100,100,100)
colour_H2 = (255,255,255)
colour_CH4 = (0,200,0)
colour_H2O = (	51	,181,	229)



class Molecule:
    def __init__(self, MASS, type, VELOCITY_X, VELOCITY_Y, POSITION_X, POSITION_Y, RADIUS, COLOUR):
        self.mass = MASS
        self.type = type
        self.v_x = VELOCITY_X
        self.v_y = VELOCITY_Y
        self.p_x = POSITION_X
        self.p_y = POSITION_Y
        self.v = math.sqrt(VELOCITY_X**2 + VELOCITY_Y**2)
        self.E_k = (1/2)*self.mass*(self.v)**2
        self.color = COLOUR
        self.radius = RADIUS
    def changeType(self, TYPE):
        self.type = TYPE
        if TYPE == "CO2":
            self.mass = mass_CO2
            self.color = colour_CO2
            self.radius = radius_CO2
        if TYPE == "H2":
            self.mass = mass_H2
            self.color = colour_H2
            self.radius = radius_H2
        if TYPE == "CH4":
            self.mass = mass_CH4
            self.color = colour_CH4
            self.radius = radius_CH4
        if TYPE == "H2O":
            self.mass = mass_H2O
            self.color = colour_H2O
            self.radius = radius_H2O

    def update(self, moleculesList):
        self.v = math.sqrt(self.v_x**2 + self.v_y**2) # update abs velocity
        self.E_k = (1 / 2) * self.mass * (self.v) ** 2 # update kinetic energy velocity

        # Check collision with molecules
        #https://en.wikipedia.org/wiki/Elastic_collision
        for otherMolecule in moleculesList:
            global n_CH4
            global n_H2
            global n_H2O
            global n_CO2
            if otherMolecule is not self:
                distanceX, distanceY = (otherMolecule.p_x - self.p_x), (otherMolecule.p_y - self.p_y)
                distance = math.sqrt(distanceX * distanceX + distanceY * distanceY)

                # check collision between two molecules
                if distance < self.radius:  # radius=3,(2radius)**2=36

                    otherMolecule.v = math.sqrt(otherMolecule.v_x ** 2 + otherMolecule.v_y ** 2)  # update abs velocity

                    phi = math.atan2(distanceY, distanceX)  # math.acos((distanceX/distance))

                    theta_self = math.atan2(self.v_y, self.v_x)
                    theta_otherMolecule = math.atan2(otherMolecule.v_y, otherMolecule.v_x)

                    self.v_y = (self.v * math.cos(theta_self - phi) *(self.mass - otherMolecule.mass) + 2*(otherMolecule.mass*otherMolecule.v*math.cos(theta_otherMolecule - phi))) / (otherMolecule.mass + self.mass) * math.cos(phi) + self.v * math.sin(theta_self - phi)*math.cos(phi*math.pi/2)
                    self.v_x = (self.v * math.cos(theta_self - phi) *(self.mass - otherMolecule.mass) + 2*(otherMolecule.mass*otherMolecule.v*math.cos(theta_otherMolecule - phi))) / (otherMolecule.mass + self.mass) * math.sin(phi) + self.v * math.sin(theta_self - phi)*math.sin(phi*math.pi/2)

                    otherMolecule.v_y = (otherMolecule.v_y *(otherMolecule.mass - self.mass) - 2*(self.mass*self.v*math.cos(theta_self - phi))) / (self.mass + otherMolecule.mass) * math.cos(phi) + otherMolecule.v * math.sin(theta_otherMolecule - phi)*math.cos(phi*math.pi/2)
                    otherMolecule.v_x = (otherMolecule.v_x *(otherMolecule.mass - self.mass) - 2*(self.mass*self.v*math.cos(theta_self - phi))) / (self.mass + otherMolecule.mass) * math.sin(phi) + otherMolecule.v * math.sin(theta_otherMolecule - phi)*math.sin(phi*math.pi/2)
                    # check if H2+H2 collision:
                    # if (self.type == "H2" and otherMolecule.type == "H2") and (self.E_k+otherMolecule.E_k > Ea_H2_H2):
                    #     self.v_x = otherMolecule.v_x
                    #     self.v_y = otherMolecule.v_y
                    # check if H2+CO2 collision
                    if ((self.type == "CO2" and otherMolecule.type == "H2") or (self.type == "H2" and otherMolecule.type == "CO")):
                        createMolecule("CH4", self.p_x, self.p_y)
                        createMolecule("H2O", otherMolecule.p_x, otherMolecule.p_y)

                        if self in moleculesList:
                            moleculesList.remove(self)
                        moleculesList.remove(otherMolecule)
                        # accounting for gained and lost molecules
                        n_CH4 += 1
                        n_H2O += 1

                        n_CO2 -= 1
                        n_H2 -= 1
                    elif (((self.type == "H2O" and otherMolecule.type == "CH4") or (self.type == "CH4" and otherMolecule.type == "H2O")) and (self.E_k + otherMolecule.E_k > Ea_H2O_CO2)):
                        createMolecule("H2", self.p_x, self.p_y)
                        createMolecule("CO2", otherMolecule.p_x, otherMolecule.p_y)
                        if self in moleculesList:
                            moleculesList.remove(self)
                        moleculesList.remove(otherMolecule)
                        # accounting for gained and lost molecules
                        n_H2 += 1
                        n_CO2 += 1

                        n_H2O -= 1
                        n_CH4 -= 1
                    self.p_y += 1*self.v_y
                    self.p_x += 1*self.v_x
                    otherMolecule.p_x += 1*otherMolecule.v_x
                    otherMolecule.p_y += 1*otherMolecule.v_y

        # Check collision with walls
        if self.p_y > HEIGHT:
            self.v_y *= -1
            self.p_y = HEIGHT

        if self.p_y < 0:
            self.v_y *= -1
            self.p_y = 0

        if self.p_x > WIDTH:
            self.v_x *= -1
            self.p_x = WIDTH

        if self.p_x < 0:
            self.v_x *= -1
            self.p_x = 0
        self.p_x += self.v_x
        self.p_y += self.v_y


def createMolecule(TYPE, P_X = -1, P_Y = -1):
    if TYPE == "CO2":
        colour = colour_CO2
        radius = radius_CO2
        mass = mass_CO2
    if TYPE == "CH4":
        colour = colour_CH4
        radius = radius_CH4
        mass = mass_CH4
    if TYPE == "H2O":
        colour = colour_H2O
        radius = radius_H2O
        mass = mass_H2O
    if TYPE == "H2":
        colour = colour_H2
        radius = radius_H2
        mass = mass_H2

    v_x = SPEED * random.random() * random.choice((-1, 1))
    v_y = SPEED * random.random() * random.choice((-1, 1))

    if (P_X == -1 and P_Y == -1):
        p_x = random.randint(0, WIDTH)
        p_y = random.randint(0, HEIGHT)
    else:
        p_x = P_X
        p_y = P_Y
    molecules_All.append(Molecule(mass, TYPE, v_x, v_y, p_x, p_y, radius, colour))

for i in range(int(n_CO2)):
    createMolecule("CO2")
for i in range(int(n_H2)):
    createMolecule("H2")
for i in range(int(n_CH4)):
    createMolecule("CH4")
for i in range(int(n_H2O)):
    createMolecule("H2O")


clock = pygame.time.Clock()


def animate(i):
    global logging_CH4
    global logging_H2O
    global logging_CO2
    global logging_H2
    global logging_time

    plt.cla()
    l = 4
    plt.plot(logging_time, logging_CH4, label='[CH4]', linewidth = l, color = (colour_CH4[0]/255, colour_CH4[1]/255, colour_CH4[2]/255))
    plt.plot(logging_time, logging_H2O, label='[H2O]', linewidth = l, color = (colour_H2O[0]/255, colour_H2O[1]/255, colour_H2O[2]/255))
    plt.plot(logging_time, logging_H2, label='[H2]', linewidth = l, color = (colour_H2[0]/255, colour_H2[1]/255, colour_H2[2]/255))
    plt.plot(logging_time, logging_CO2, label='[CO2]', linewidth = l, color = (colour_CO2[0]/255, colour_CO2[1]/255, colour_CO2[2]/255))

    plt.legend(loc='upper right', fontsize=14)
    plt.tight_layout()

logging_CH4 =[]
logging_H2O =[]
logging_H2 = []
logging_CO2 = []
logging_time = []
time = 0

# draw molecules
fig = plt.figure()
axes = fig.add_subplot(111)
line, = plt.plot(logging_time, logging_H2O)

ani = FuncAnimation(plt.gcf(), animate, interval=1000)
plt.tight_layout()
plt.xlabel("Time (s)")
plt.ylabel("Concentration (Mol/L)")
pygame.time.wait(5000)

while True:

    for event in pygame.event.get():
        if event.type == pygame.KEYDOWN:
            # As long as an arrow key is held down, the respective speed is set to 3 (or minus 3)
            if event.key == pygame.K_x:
                for z in range(add):
                    createMolecule("H2", random.randint(0,WIDTH), random.randint(0,HEIGHT))
                    n_H2 += 1
            if event.key == pygame.K_RIGHT:
                showGraph = not showGraph
            if event.key == pygame.K_UP:
                tickrate *= 2
            if event.key == pygame.K_DOWN:
                tickrate /= 2
    screen.fill((255, 157, 111))

    for molecule in molecules_All:
        molecule.update(molecules_All)
        pygame.draw.circle(screen, molecule.color, (int(molecule.p_x), int(molecule.p_y)), int(molecule.radius))

    clock.tick(tickrate)
    time += clock.get_rawtime()/1000
    pygame.display.flip()

    logging_CH4.append(n_CH4/V)
    logging_H2O.append(n_H2O/V)
    logging_H2.append(n_H2/V)
    logging_CO2.append(n_CO2/V)
    logging_time.append(time)

    plt.draw()
    line.set_data(logging_time, logging_H2O)
    axes.autoscale_view(True,True,True)
    axes.relim()
    plt.draw()

    if showGraph:
        plt.pause(0.000001)


