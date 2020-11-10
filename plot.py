import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.ndimage.filters import gaussian_filter1d
from matplotlib import animation

#  plt.rcParams['animation.ffmpeg_path'] = ''

# y = np.fromfile('out/A')
# y2 = np.fromfile('out/A')
# y3 = np.fromfile('out/A')
SAVE_NB = 653

y = np.fromfile('res/A_A')
y2 = np.fromfile('res/B_A')
y3 = np.fromfile('res/B_A')
y = np.reshape(y, (SAVE_NB, y.shape[0] // SAVE_NB))
y2 = np.reshape(y2, (SAVE_NB, y2.shape[0] // SAVE_NB))
y3 = np.reshape(y3, (SAVE_NB, y3.shape[0] // SAVE_NB))
#  y = np.fromfile('res/A')
#  y2 = np.fromfile('res/A')
#  y3 = np.fromfile('res/A')
#  from sys import argv
#  SAVE_NB = int(argv[1])
#  # SAVE_NB = 87
#  y = np.reshape(y, (SAVE_NB, y.shape[0] // SAVE_NB))
#  # SAVE_NB = 173
#  y2 = np.reshape(y2, (SAVE_NB, y2.shape[0] // SAVE_NB))
#  # SAVE_NB = 345
#  y3 = np.reshape(y3, (SAVE_NB, y3.shape[0]//SAVE_NB))

# y = y[150:, :]
# y2 = y2[150:, :]
# y3 = y3[150:, :]


print(y.shape)
print(y2.shape)
print(y3.shape)

# print(np.max(np.abs(y[241, :]-y2[241, :])))
# print(np.max(np.abs(y[241, :])))

# print(y[241,:])

# print((y-y2).sum().sum())
# print((y-y3).sum().sum())
# print((y2-y3).sum().sum())

# print(y[0])
from scipy import interpolate

#  y = y / y[0]
#  y2 = y2 / y2[0]
#  y3 = y3 / y3[0]

tot = len(y[0]) + len(y2[0]) + len(y3[0])

FRAME_NB = y.shape[0]

y = y - y[0]
y2 = y2 - y2[0]
y3 = y3 - y3[0]

fig, ax = plt.subplots(figsize=(16, 9))
x = np.linspace(0, 3, len(y[0]))
x2 = np.linspace(3, 6, len(y2[0]))
x3 = np.linspace(3, 6, len(y3[0]))
line, = ax.plot(x, np.sin(x / 100), 'b', label='A')
line2, = ax.plot(x2, np.sin(x2 / 100), 'r', label='B')
line3, = ax.plot(x3, np.sin(x3 / 100), 'g')
plt.legend()
time_text = ax.text(0.05, 0.95, '', horizontalalignment='left',
                    verticalalignment='top', transform=ax.transAxes)

print(np.min(np.min(y)))
print(np.max(np.max(y)))
ax.set_ylim([np.min([np.min(np.min(y)), np.min(np.min(y2)), np.min(np.min(y3))]),
             np.max([np.max(np.max(y)), np.max(np.max(y2)), np.max(np.max(y3))])])
# ax.set_ylim([0.0, 0.06])

# ypsmooth = gaussian_filter1d(y[i], sigma=3)


def animate(i):
    line.set_ydata(y[i])  # update the data.
    line2.set_ydata(y2[i])  # update the data.
    line3.set_ydata(y3[i])  # update the data.
    time_text.set_text('frame ' + str(i) + '/' + str(FRAME_NB))

    # return line, time_text
    #  return line,  time_text
    return line, line2, line3, time_text


ani = FuncAnimation(fig, animate, frames=len(y[:, 0]), interval=1, blit=True)

plt.show()

mywriter = animation.FFMpegWriter(fps=60, codec='libx264', extra_args=[
                                 '-pix_fmt', 'yuv420p', '-profile:v', 'high', '-tune', 'animation', '-crf', '18'])

ani.save('seawave_1d_ani.mp4', writer=mywriter)
