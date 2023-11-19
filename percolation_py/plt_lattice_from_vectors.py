import matplotlib.pyplot as plt
import numpy as np

vectors = [np.array([1, 0, 1]), np.array([1, -1, -1]), np.array([1, 1, 0]), np.array([1, 1, -1])]

box = [3, 3, 3]
x_list = []
y_list = []
z_list = []

for x in range(box[0]):
    for y in range(box[1]):
        for z in range(box[2]):
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for x in range(box[0]):
    for y in range(box[1]):
        for z in range(box[2]):
            for v in range(len(vectors)):
                if 0 <= x + vectors[v][0] < box[0] and 0 <= y + vectors[v][1] < box[1] and 0 <= z + vectors[v][2] < box[2]:
                    ax.plot([x, x + vectors[v][0]], [y, y + vectors[v][1]], [z, z + vectors[v][2]], '-k')
                if 0 <= x - vectors[v][0] < box[0] and 0 <= y - vectors[v][1] < box[1] and 0 <= z - vectors[v][2] < box[2]:
                    ax.plot([x, x - vectors[v][0]], [y, y - vectors[v][1]], [z, z - vectors[v][2]], '-k')
col = 'rgb'
for k in range(len(x_list)):
    ax.plot([x_list[k]], [y_list[k]], [z_list[k]], 'o'+col[z_list[k]]) # plot nodes
ax.set(xticklabels=[], yticklabels=[], zticklabels=[])
ax.grid(False)
ax.tick_params(left=False, bottom=False)
ax.view_init(elev=20, azim=5)
plt.axis('off')
fig.set_size_inches(4, 4)
fig.savefig('3D_optim_lattice.pdf')
plt.show()

vectors = [np.array([5, 7]), np.array([7, 4]), np.array([0, 3]), np.array([7, -6]), np.array([7, -5])]
fig = plt.figure()
ax = fig.add_subplot()
ax.plot([i % 15 - 7 for i in range(15*15)], [np.floor(i/15) - 7 for i in range(15*15)], 'ok')
for k in range(len(vectors)):
    ax.plot([vectors[k][0], -vectors[k][0]], [vectors[k][1], -vectors[k][1]], '-k')
plt.axis('off')
fig.set_size_inches(3.3, 3.3)
ax.set_aspect('equal', adjustable='box')
fig.savefig('2D_optim_lattice_old.pdf')
plt.show()
