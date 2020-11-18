import getopt
import numpy as np
import matplotlib.pyplot as plt
import math
import json
import scipy.interpolate as interpolate

plt.style.use(['science', 'grid'])


class Simulation:
    def __init__(self, path, quantities, openbf=False):
        if not path.endswith('/'):
            path += '/'
        with open(path + 'metadata.json') as f:
            data = json.load(f)
        vessels_infos = data['vessels_info']
        self.save_nb = data['save_nb']
        self.time_simul = data['total_time']
        self.path = path
        self.x = np.linspace(0.0, self.time_simul, self.save_nb, endpoint=True)
        self.dx = data['dx']
        self.vessels = []

        for name,length in vessels_infos:
            self.vessels.append(Vessel(path, name, length, self.save_nb, quantities, self.x, openbf=openbf))

        self.dx = self.vessels[0].length / self.vessels[0].data[quantities[0]].shape[1]

    def post_process(self, quantity_to_process, *funcs, **kwargs):
        for v in self.vessels:
            v.post_process(quantity_to_process, *funcs, dx=self.dx, **kwargs)

    def plot_quantity(self, quantity_to_plot, name, xlabel='Time ($s$)', ylabel='', title=None):
        fig, ax = plt.subplots(figsize=(16,9))

        for v in self.vessels:
            for quantity,y in v.data_processed.items():
                if quantity != quantity_to_plot:
                    continue
                label = '{} ({})'.format(v.name, quantity)
                ax.plot(self.x, y, label=label, linewidth=2)

        ax.legend(fontsize=15)
        plt.xlabel(xlabel, fontsize=15)
        plt.ylabel(ylabel, fontsize=15)
        if title != None:
            plt.title(title, fontsize=18)

        fig.savefig(self.path + name, bbox_inches='tight')

    def plot_quantity_compare(self, sim, ind, quantity_to_plot, name, xlabel='Time ($s$)', ylabel='', title=None):
        fig, ax = plt.subplots(figsize=(16,9))

        #  for v in self.vessels:
        v = self.vessels[ind]
        for quantity,y in v.data_processed.items():
            if quantity != quantity_to_plot:
                continue
            label = '{} ({}) 1'.format(v.name, quantity)
            ax.plot(self.x, y, label=label, linewidth=2)

        v = sim.vessels[ind]
        for quantity,y in v.data_processed.items():
            if quantity != quantity_to_plot:
                continue
            label = '{} ({}) openbf'.format(v.name, quantity)
            ax.plot(sim.x, y, label=label, linewidth=2)

        ax.legend(fontsize=15)
        plt.xlabel(xlabel, fontsize=15)
        plt.ylabel(ylabel, fontsize=15)
        if title != None:
            plt.title(title, fontsize=18)

        fig.savefig(self.path + name, bbox_inches='tight')

class Vessel:
    def __init__(self, path, name, length, save_nb, quantities, x, openbf=False):
        self.name = name
        self.length = length
        self.data = {}
        self.data_processed = {}

        for quantity in quantities:
            if openbf:
                y = np.loadtxt(path + name + '_' + quantity)
                oldx = y[:,0]
                y = y[:,1:]
                f = interpolate.interp1d(oldx, y, axis=0, kind='cubic')
                y = f(x)
            else:
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

results_path = 'res'
sim = Simulation(results_path, ['A', 'u'], openbf=False)

#  results_path = 'res_openbf'
#  sim_openbf = Simulation(results_path, ['A', 'u'], openbf=True)

#  results_path = 'res2'
#  sim_openbf = Simulation(results_path, ['A', 'u'], openbf=False)

sim.post_process('A', get_data_at_position, pos=0.2)
#  sim_openbf.post_process('A', get_data_at_position, pos=0.2)
sim.plot_quantity('A', 'u_diff6.png', ylabel=r'$A - A_0$', title='A, x = 6m')

#  sim.plot_quantity_compare(sim_openbf, 0, 'A', 'u_diff6.png', ylabel=r'$A - A_0$', title='A, x = 6m')
#
#  sim.post_process('A', get_data_at_position, pos=2.0)
#  sim_openbf.post_process('A', get_data_at_position, pos=2.0)
#  sim.plot_quantity_compare(sim_openbf, 0, 'A', 'u_diff.png', ylabel=r'$A - A_0$', title='A, x = 2m')
#
#  sim.post_process('A', get_data_at_position, pos=0.0)
#  sim_openbf.post_process('A', get_data_at_position, pos=0.0)
#  sim.plot_quantity_compare(sim_openbf, 0, 'A', 'A_diff.png', ylabel=r'$A - A_0$', title='A, x = 0m')

#  sim.post_process('A', get_data_at_position, relative_change, pos=4.0)
#  sim_openbf.post_process('A', get_data_at_position, relative_change, pos=4.0)
#  sim.plot_quantity_compare(sim_openbf, 'A', 'A_relative.png', ylabel=r'$\frac{A - A_0}{A_0}$', title='A, x = 4m')
#
#  sim.post_process('u', get_data_at_position, pos=2.0)
#  sim_openbf.post_process('u', get_data_at_position, pos=2.0)
#  sim.plot_quantity_compare(sim_openbf, 'u', 'u.png', ylabel=r'$u$', title='u, x = 2m')

