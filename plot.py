import getopt
import numpy as np
import matplotlib.pyplot as plt
import math
import json

plt.style.use(['science', 'grid'])

results_path = 'res'

class Simulation:
    def __init__(self, path, quantities):
        if not path.endswith('/'):
            path += '/'
        with open(path + 'metadata.json') as f:
            data = json.load(f)
        vessels_infos = data['vessels_info']
        self.dx = data['dx']
        self.save_nb = data['save_nb']
        self.time_simul = data['total_time']
        self.path = path
        self.x = np.linspace(0.0, self.time_simul, self.save_nb, endpoint=True)
        self.vessels = []

        for name,length in vessels_infos:
            self.vessels.append(Vessel(path, name, length, self.save_nb, quantities))

    def post_process(self, quantity_to_process, *funcs, **kwargs):
        for v in self.vessels:
            v.post_process(quantity_to_process, *funcs, dx=self.dx, **kwargs)

    def show_quantity(self, quantity_to_plot, name, xlabel='Time ($s$)', ylabel='', title=None):
        fig, ax = plt.subplots(figsize=(16,9))

        self.plot_quantity(ax, quantity_to_plot)

        ax.legend(fontsize=15)
        plt.xlabel(xlabel, fontsize=15)
        plt.ylabel(ylabel, fontsize=15)
        if title != None:
            plt.title(title, fontsize=18)

        fig.savefig(self.path + name, bbox_inches='tight')

    def plot_quantity(self, ax, quantity_to_plot):
        for v in self.vessels:
            for quantity,y in v.data_processed.items():
                if quantity != quantity_to_plot:
                    continue
                label = '{} ({})'.format(v.name, quantity)
                ax.plot(self.x, y, label=label, linewidth=2)
    
class Vessel:
    def __init__(self, path, name, length, save_nb, quantities):
        self.name = name
        self.length = length
        self.data = {}
        self.data_processed = {}
        for quantity in quantities:
            y = np.fromfile(path + name + '_' + quantity)
            y = np.reshape(y, (save_nb, y.shape[0] // save_nb))
            self.data[quantity] = y

    def post_process(self, quantity_to_process, *funcs, **kwargs):
        for quantity,y in self.data.items():
            if quantity != quantity_to_process:
                continue
            tmp = y
            for func in funcs:
                tmp = func(tmp, **kwargs)
            self.data_processed[quantity] = tmp



def get_data_at_position(data, **kwargs):
    dx = kwargs.get('dx', None)
    pos = kwargs.get('pos', None)
    x = pos / dx
    xa = math.floor(x)
    xb = xa + 1
    ya = data[:, xa]
    yb = data[:, xb]
    slope = yb - ya
    yc = ya + (x - xa) * slope
    return yc

def difference(data, **kwargs):
    return data - data[0]

def relative_change(data, **kwargs):
    return difference(data) / data[0]



sim = Simulation(results_path, ['A', 'u'])

sim.post_process('A', get_data_at_position, difference, pos=2.0)
sim.show_quantity('A', 'A_diff.png', ylabel=r'$A - A_0$', title='A, x = 2m')

sim.post_process('A', get_data_at_position, relative_change, pos=4.0)
sim.show_quantity('A', 'A_relative.png', ylabel=r'$\frac{A - A_0}{A_0}$', title='A, x = 4m')

sim.post_process('u', get_data_at_position, pos=2.0)
sim.show_quantity('u', 'u.png', ylabel=r'$u$', title='u, x = 2m')

