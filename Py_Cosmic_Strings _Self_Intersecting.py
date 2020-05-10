""""
Simulating network of cosmic strings with fixed ends on cubic surface.
//
Uses algorithm adapted from Smith & Vilenkin (Physical Review D, Vol 36, Number 4)
//
String class has a random walk function which it uses to initialize itself, then
dynamics are handled by finite difference method for wave equation.
//
author: Christian Capanelli
christiancapanelli@gmail.com
"""

# Imports
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from matplotlib import rc

# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Computer Modern Roman']})
# rc('text', usetex=True)

# Box parameters
N = 32 # number length of box
delta = 1 # step-size
boxL = N * delta  # enclosing box length
rep = 8  # How many times random walk steps are repeated

self_intersect = False
save_anim = False


# -- Define CosmicString Class -- #
class CosmicString:

    def __init__(self, seed, time_param, repeat_param):
        self.repeat = repeat_param
        self.worldsheet = np.empty((1, 1))  # Initialize worldsheet to be overwritten later
        self.T = time_param  # accepts end time to integrate to
        self.seed = seed
        np.random.seed(self.seed)

        # Set origin of walk inside given box:
        x0 = int(np.random.uniform(-boxL / 2, boxL / 2))
        y0 = int(np.random.uniform(-boxL / 2, boxL / 2))
        z0 = int(np.random.uniform(-boxL / 2, boxL / 2))
        origin = [x0, y0, z0]

        self.moves = delta * np.array([[0, 0, 1], [0, 1, 0],
                                       [1, 0, 0], [-1, 0, 0],
                                       [0, 0, -1], [0, -1, 0]])  # the possible moves for walk
        self.V_cusp = self.moves / delta  # Defining possible initial velocities
        self.V_diag = (1 / 2) * np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1],
                                          [-1, -1, 0], [0, -1, -1], [-1, 0, -1],
                                          [1, -1, 0], [-1, 1, 0], [0, 1, -1],
                                          [0, -1, 1], [1, 0, -1], [-1, 0, 1]])
        self.length = 1
        self.loops = []
        self.walk(origin)

    def walk(self, origin):
        x = 0
        y = 0
        z = 0

        i = 1
        coordinates = [origin]
        # Walks until box edge is hit
        while np.abs(x) < boxL / 2 and np.abs(y) < boxL / 2 and np.abs(z) < boxL / 2:
            step = self.moves[np.random.randint(0, 6)]

            for k in range(self.repeat):
                x = coordinates[i - 1][0] + step[0]
                y = coordinates[i - 1][1] + step[1]
                z = coordinates[i - 1][2] + step[2]

                coord = [x, y, z]
                coordinates.append(coord)

                i += 1
        self.length = i
        self.worldsheet = np.empty((int(self.T / delta) + 1, self.length, 3))
        self.worldsheet[0, :, :] = coordinates

    def dynamics(self):
        t = 0  # actual time
        j = 0  # time index
        while t < self.T:  # loop through time
            coordinates = np.empty((self.length, 3))
            for i in range(0, self.length):  # loop along spacelike parameter at single timeslice

                # Fix end points:
                if i == 0 or i == self.length - 1 or np.isnan(self.worldsheet[j, i, 0]) or np.isnan(
                        self.worldsheet[j, i + 1, 0]) or np.isnan(self.worldsheet[j, i - 1, 0]):
                    X_next = self.worldsheet[j, i, :]
                    coordinates[i, :] = X_next
                else:
                    if j == 0:  # instruction for time zero
                        X_current = self.worldsheet[j, i, :]
                        X_back = self.worldsheet[j, i - 1, :]
                        X_forward = self.worldsheet[j, i + 1, :]
                        if np.abs(np.linalg.norm(X_forward - X_back)) <= 1e-3:
                            X_next = X_current + delta * self.V_cusp[np.random.randint(0, 6)]
                        elif np.abs(np.linalg.norm(X_forward - X_back) - 1) <= 1e-3:
                            X_next = X_current + delta * (self.V_diag[np.random.randint(0, 12)])
                        else:
                            X_next = X_current
                    else:
                        k = i + 1
                        while np.isnan(self.worldsheet[j, k, 0]):
                            k += 1
                            # if k == self.length :
                            #     break
                        X_back = self.worldsheet[j, i - 1, :]
                        X_forward = self.worldsheet[j, k, :]
                        X_prev = self.worldsheet[j - 1, i, :]
                        X_next = X_forward + X_back - X_prev

                    coordinates[i, :] = X_next

                    # Check for self-intersection:
                    if self_intersect:
                        for m in range(i + 1, self.length):
                            rtol = 1e-05
                            atol = 1e-08
                            if np.allclose(self.worldsheet[j, i, :], self.worldsheet[j, m, :], rtol, atol):
                                inherited_time = t
                                inherited_time_index = j
                                time_param = self.T
                                inherited_worldsheet = self.worldsheet[:, i:m + 1, :]
                                self.loops.append(CosmicStringLoop(inherited_time, inherited_time_index, time_param,
                                                                   inherited_worldsheet))
                                coordinates[i:m + 1, :] = np.array([np.nan, np.nan, np.nan])

            self.worldsheet[j + 1, :, :] = coordinates[:, :]
            t += delta
            j += 1


class CosmicStringLoop:
    def __init__(self, inherited_time, inherited_time_index, time_param, inherited_worldsheet):
        self.time_born = inherited_time
        self.birth_index = inherited_time_index
        self.worldsheet = inherited_worldsheet  # Takes sliced worldsheet of parent string
        self.T = time_param  # accepts end time to integrate to
        self.length = np.shape(self.worldsheet)[1]

    def dynamics(self):
        while self.time_born < self.T:  # loop through time
            coordinates = np.empty((self.length, 3))
            for i in range(0, self.length):  # loop along spacelike parameter at single timeslice

                X_back = self.worldsheet[self.birth_index, i - 1, :]
                X_forward = self.worldsheet[self.birth_index, i + 1, :]
                X_prev = self.worldsheet[self.birth_index - 1, i, :]
                X_next = X_forward + X_back - X_prev

                coordinates[i, :] = X_next  # String coordinates at single time slice
            self.worldsheet[self.birth_index + 1, :, :] = coordinates[:, :]
            self.time_born += delta
            self.birth_index += 1


# -- Initialize some cosmic strings -- #

number_strings = 5
number_frames = 50
end_time = number_frames * delta

strings = []

# Pre-calculate dynamics for given end time
for n in range(number_strings):
    strings.append(CosmicString(n, end_time, rep))
    strings[n].dynamics()

## -- VISUALS -- #
framerate = 12

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111, projection='3d')


# -- ANIMATION BLOCK -- #
def animate(i):
    ax.clear()
    plt.title(
        r"$V= ($" + str(N) + "$\delta)^3$; Number of Strings = "
        + str(number_strings), fontdict={'fontsize': 11})

    ax.view_init(45, -45)

    ax.set_xlim3d([-boxL / 2, boxL / 2])
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['linewidth'] = 0.1
    ax.set_xticklabels([])

    ax.set_ylim3d([-boxL / 2, boxL / 2])
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis._axinfo["grid"]['linewidth'] = 0.1
    ax.set_yticklabels([])

    ax.set_zlim3d([-boxL / 2, boxL / 2])
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis._axinfo["grid"]['linewidth'] = 0.1
    ax.set_zticklabels([])

    ax.text2D(0.8, 0.1, str(i) + r"$\delta$", transform=ax.transAxes)

    for n in range(number_strings):
        string_plot(strings, i, n)
        for k in range(np.shape(strings[n].loops)[0]):
            ax.plot(strings[n].loops[k].worldsheet[i, :, 0],
                    strings[n].loops[k].worldsheet[i, :, 1],
                    strings[n].loops[k].worldsheet[i, :, 2],
                    linewidth=0.3, color='k')


def string_plot(strings_arg, i, j):
    return ax.plot(strings_arg[j].worldsheet[i, :, 0],
                   strings_arg[j].worldsheet[i, :, 1],
                   strings_arg[j].worldsheet[i, :, 2],
                   linewidth=0.3, color='k')


anim = animation.FuncAnimation(fig, animate, frames=number_frames)  # Call the animator

if save_anim:
    anim.save('gauged_string_intersection_animation_' + str(number_strings) + '_strings_N' + str(N) + '.mp4',
              fps=framerate, extra_args=[
            '-vcodec', 'libx264'])

###################################################################

plt.show()
