import os
from conditions import initial_conditions
from animate import animate
import multiprocessing as mp

os.system("g++ -std=c++17 vector3d.cpp body.cpp simulation.cpp -o sim")

T = 100000
Tsave = 1000
dt = 1e-3

N = [1]
L = [1]

def simulate(n,l):
    initial_conditions(n , l)

    output_file = f"{n}particles_{l}m"
    input_file = f"random_particles_{n}particles_{l}m.csv"


    os.system("(cd ..; chmod 700 -R assessment4; cd assessment4)")
    os.system(f"./sim cache/{input_file} cache/{output_file}.csv {dt} {T} {Tsave} {l}")
    animate(output_file, l)
    # os.system(f"python3 animate.py cache/{output_file}.csv animations/{output_file}.gif {l}")


if __name__ == '__main__':
    for n in N:
        for l in L:
            test = mp.Process(target=simulate, args=(n,l))
            test.start()


# n = int(input("Enter number of particles: "))
# L = float(input("Size of the box (meters): "))
#
# initial_conditions(n,L)
#
# T = input("Enter number of steps for simulation: ")
# Tsave = input("Enter the save interval (lowest is 1 = every step)(lower = more ram but better animations): ")
# dt = input("Enter dt (smaller step = better resolution = longer runtime): ")

# input_file = "cache/random_particles.csv"
# output_file = f"{n}particles_{L}m"
#
# os.system("(cd ..; chmod 700 -R assessment4; cd assessment4)")
# os.system("g++ -std=c++17 vector3d.cpp body.cpp simulation.cpp -o sim")
# os.system(f"./sim {input_file} cache/{output_file}.csv {dt} {T} {Tsave} {L}")
# os.system(f"python3 animate.py cache/{output_file}.csv animations/{output_file}.gif {L}")
