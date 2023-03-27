import csv
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation
import numpy as np

L = 1
L /= 2

savefile = "test_case2.csv"

file = open(savefile)
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


# Since plotting a single graph
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1,projection="3d")


def animation_function(i):
    ax.cla()
    ax.set_xlim([-L, L])
    ax.set_ylim([-L, L])
    ax.set_zlim([-L, L])
    ax.set_title('Particle Motion')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
    for n in range(N):
        frame = ax.scatter(pos[i, n, 0], pos[i, n, 1], pos[i, n, 2])
    return frame
#
# def generate_points():
# 	while
anim = FuncAnimation(fig,
                    func = animation_function,
                    frames = NumSave,
                    interval = 0.1,
                    repeat = True)

writergif = animation.PillowWriter(fps = 30)
anim.save("motion_ani.gif", writer=writergif)
