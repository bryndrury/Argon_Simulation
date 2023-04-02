import numpy as np

def initial_conditions(n, L):
    # number of particles
    # size of the box (and frame)

    index = np.linspace(1,n,n)
    # mass of argon
    mass = np.ones(n*3)
    mass *= 6.6335209e-26

    pos = np.random.uniform(-L/2, L/2, n*3)
    vel = np.random.uniform(-L/4, L/4, n*3)
    # vel = np.zeros(n*3)

    # pos = [4.35e-7,0,0,-4.35e-7,0,0]
    # vel = np.zeros(6)
    file = open(f"cache/random_particles_{n}particles_{L}m.csv", "w+")

    for i in range(0,n):
        string = ["Particle "+str(int(index[i]))+"\n",str(mass[i])+"\n",str(pos[3*i])+", "+str(pos[3*i+1])+", "+str(pos[3*i+2])+"\n",str(vel[3*i])+", "+str(vel[3*i+1])+", "+str(vel[3*i+2])+"\n"]
        file.writelines(string)
    file.close()
