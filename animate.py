import csv
import os
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation
import numpy as np
import sys

# n = len(sys.argv)
#
# print("/nRunning Animation.")
# print(f"Input file: {sys.argv[1]}")
# print(f"Output file: {sys.argv[2]}")
# print(f"Size of Box: {sys.argv[3]}")
#
# if n != 4:
#     print("Invalid input.")
#     print("animate <input file.csv> <output file.gif> <L>")
#     exit()
# savefile = str(sys.argv[1])
#
# L = float(sys.argv[3])
def animate(input_file, L):
    savefile = f"cache/{input_file}.csv"
    output_file = f"animations/{input_file}.gif"

    try:
        file = open(savefile)
    except:
        print("Failed to open file.")
        # exit()
    csvreader = csv.reader(file)
    rows = []
    for row in csvreader:
        rows.append(row)
    file.close()

    dt = float(rows[0][0])
    T = int(rows[0][1])
    Tsave = int(rows[0][2])
    N = int(rows[0][3])
    NumSave = int(np.ceil((T+1)/Tsave))

    names = rows[1][:]

    time = np.zeros([NumSave,1])
    pos = np.zeros([NumSave,N,3])

    data = np.zeros(3*N)
    for i in range(NumSave):
        start = i * (N+1) + 2
        time[i] = float(rows[start][0])

        for j in range(N):
            for k in range(3):
                pos[i, j, k] = float(rows[start + j + 1][k])






    # Setting limits for x and y axis
    print("Animating... ")

    # Since plotting a single graph
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1,projection="3d")


    def animation_function(i):
        ax.cla()
        ax.set_xlim([-L, L])
        ax.set_ylim([-L, L])
        ax.set_zlim([-L, L])
        ax.set_title(f'Particle Motion, {N} particles')
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_zlabel('$z$')
        for n in range(N):
            frame = ax.scatter(pos[i, n, 0], pos[i, n, 1], pos[i, n, 2])
        os.system("clear")
        print(f"{round(100*i/NumSave, 1)}%")
        return frame


    anim = FuncAnimation(fig,
                        func = animation_function,
                        frames = NumSave,
                        interval = 0.1,
                        repeat = True)

    writergif = animation.PillowWriter(fps = 30)
    os.system("clear")
    print("100%")
    anim.save(f"{output_file}", writer=writergif, dpi=200)
    # anim.save(f"{sys.argv[2]}", writer=writergif, dpi=200)
