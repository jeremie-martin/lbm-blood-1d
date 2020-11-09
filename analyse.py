import numpy as np
import matplotlib as mpl
# mpl.use("pgf")
import matplotlib.pyplot as plt
# plt.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })

# import matplotlib.pyplot as plt
import sys

y = np.loadtxt('res/A')
y2 = np.loadtxt('res/B')
y3 = np.loadtxt('res/C')

# y = np.fromfile('A4008')
# y2 = np.fromfile('A4004')
# y3 = np.fromfile('A4002')
# y3 = np.loadtxt('A1_A.out')

# SAVE_NB = 139
# y = np.reshape(y, (SAVE_NB, y.shape[0] // SAVE_NB))
# SAVE_NB = 277
# y2 = np.reshape(y2, (SAVE_NB, y2.shape[0] // SAVE_NB))
# SAVE_NB = 553
# y3 = np.reshape(y3, (SAVE_NB, y3.shape[0]//SAVE_NB))

#  y = np.fromfile('out/A')
#  y2 = np.fromfile('out/A')
#  y3 = np.fromfile('out/A')
#
#  # y3 = np.loadtxt('A1_A.out')
#
#  SAVE_NB = 144
#  y = np.reshape(y, (SAVE_NB, y.shape[0] // SAVE_NB))
#  SAVE_NB = 144
#  y2 = np.reshape(y2, (SAVE_NB, y2.shape[0] // SAVE_NB))
#  SAVE_NB = 144
#  y3 = np.reshape(y3, (SAVE_NB, y3.shape[0]//SAVE_NB))

# y = np.loadtxt('owo.out')

# x1 = y[:, 0]
# y = y[:, 1:]

# SAVE_NB = 5001
# SAVE_NB2 = int(sys.argv[1])
# y = np.reshape(y, (SAVE_NB, y.shape[0]//SAVE_NB))
# y2 = np.reshape(y2, (SAVE_NB2, y2.shape[0]//SAVE_NB2))

# y = np.loadtxt("A1_A.out")


# idx = (np.abs(x3 - (i*(3.3/SAVE_NB)))).argmin()

ttt = int(1*y.shape[0]/3)
ttt2 = int(1*y2.shape[0]/3)
ttt3 = int(1*y3.shape[0]/3)
(ttt, ttt2, ttt3) = (0,0,0)
print(ttt, ttt2, ttt3)
print(y.shape, y2.shape, y3.shape)
y1 = y[ttt:, int(y.shape[1]/2)-1]
y2 = y2[ttt2:, int(y2.shape[1]/2)-1]
y3 = y3[ttt3:, int(y3.shape[1]/2)-1]
#  y1 = y[:, 0]
#  y2 = y2[:, 0]
#  y3 = y3[:, 0]

fig, ax = plt.subplots(figsize=(11.69,8.27))
x1 = np.linspace(0, 3.3, len(y1), endpoint=True)
x2 = np.linspace(0, 3.3, len(y2), endpoint=True)
x3 = np.linspace(0, 3.3, len(y3), endpoint=True)

print(np.min(y1), np.max(y1))
print(np.min(y2), np.max(y2))
#  y1 = (y1 / y1[0])
#  y2 = y2 / y2[0]
#  y3 = y3 / y3[0]
print(np.min(y1), np.max(y1))
print(np.min(y2), np.max(y2))


ax.tick_params(axis='both', which='major', labelsize=13)
ax.tick_params(axis='both', which='minor', labelsize=13)

# ax.semilogy(x, yStart - yEnd)
ax.plot(x1, y1,  label='$\Delta x = 0.01$', linewidth=3)
ax.plot(x2, y2,  label='$\Delta x = 0.0075$', linewidth=3)
ax.plot(x3, y3, label='$\Delta x = 0.005$',  linewidth=3)

ax.legend(fontsize=15)

plt.title('$L = 6m$, $x = 3m$', fontsize=18)

plt.xlabel('Time ($s$)', fontsize=15)
plt.ylabel('$\\frac{A}{A_0}$', fontsize=15)

# fig.savefig("after.png", bbox_inches='tight')
# plt.savefig("output2.pdf")

# ax.axvline(x=x[startMaxIdx])
# ax.axvline(x=x[endMaxIdx], color='r')
plt.show()
print('lol')
