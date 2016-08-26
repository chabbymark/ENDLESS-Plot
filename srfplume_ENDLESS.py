#!/usr/bin/env python3

# Program: check mean surface plume for ENDLES
# 06/30/2015  Bicheng Chen  First released

### Import module ###
# module
import sys
sys.path.append("../")# use to import modules from parent directory
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib import cm
from glob import glob
import struct
import numpy as np
import lespy as lp
from numpy import pi
# environment
### END Import module ###

### User specified Variables ###
# File
path = '/data/1/bzc/LES/Coupling/LV3_fromYang/output'
ftype = 'con'

# simulation
tt = 280000
dm = lp.dmClass.domain(nx=100, ny=100, nz=151,
    dx=5, dy=5, dz=1, lx=500, ly=500, lz=150)
dm.show()
zi = 100
ustar = 0.0061

# figure
sz = dict(left=0.1, right=0.88, bottom=0.1, top=0.98,
    wspace=0.05, hspace=0.05)
level_bd = (-8, -4)
xlim = (0, 2)
ylim = (0, 1)
lw = 0.7
ars = 8

### END User specified Variables ###

### Plot Variables ###
con = dict()
#u_cell = dict()
#v_cell = dict()
### END Plot Variables ###

### Predefined Function ###
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
  new_cmap = colors.LinearSegmentedColormap.from_list(
      'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
      cmap(np.linspace(minval, maxval, n)))
  return new_cmap
### END Predefined Function ###

### Main Body ###
## Read file
fn_pre = ftype + '_tt' + str(tt).zfill(8)
fn_ls = sorted(glob(path + '/' + fn_pre + '*'))
id_con = list(set([fn_ele.split('_s')[-1].split('_r')[0] for fn_ele in fn_ls]))
fn_con = dict()
for ind in range(len(id_con)):
  fn_con[ind] = [fn_ele for fn_ele in fn_ls if ('s' + str(id_con[ind]).zfill(3)) in fn_ele]
nd = dm.Ld*dm.Ny*dm.Nz
for ind in range(len(id_con)):
  con[id_con[ind]] = dict()
  for fn in fn_con[ind]:
    row, col = fn.split('_r')[-1].split('_c')
    col = col.split('.out')[0]
    if row not in con[id_con[ind]].keys():
      con[id_con[ind]][row] = dict()
    with open(fn, 'rb') as f:
      f.read(4)
      con[id_con[ind]][row][col] =\
          np.asfarray(struct.unpack('d'*nd, f.read(8*nd)),
          dtype='float').reshape((dm.Nz, dm.Ny, dm.Ld))[0,:,0:dm.Nx]

## Process data
## Plot the data
# initialization
levels_con = np.arange(level_bd[0],level_bd[1]+0.1,0.1)
cmap_con = cm.get_cmap("jet")
#cmap_con = cm.get_cmap("RdYlBu_r")
#cmap_con = cm.get_cmap("Greys_r")
#cmap_con = truncate_colormap(cm.get_cmap("Greys_r"),0.2, 1)
ifig = 0
# plot the instaneous surface concentration
for ind in range(len(id_con)):
  ifig += 1
  fig = plt.figure(ifig, figsize=[8, 5])
  gs = gridspec.GridSpec(1,1)
  gs.update(**sz)
  ax = plt.subplot(gs[0, 0])
  ax.set_aspect('equal', adjustable='box')
  for row in con[id_con[ind]].keys():
    for col in con[id_con[ind]][row].keys():
      row_i = int(row)
      col_i = int(col)
      y, x = np.mgrid[(row_i-1)*dm.Ly:row_i*dm.Ly:dm.Dy,
          (col_i-1)*dm.Lx:col_i*dm.Lx:dm.Dx] / 1000
      con_plt = np.array(con[id_con[ind]][row][col])
      print(row, col)
      print(id_con[ind], np.amax(con_plt))
      idx, idy = np.where(con_plt < 1e-50)
      con_plt[idx, idy] = 1e-50
      caxc = ax.contourf(x, y, np.log10(con_plt),
          levels_con, cmap=cmap_con, extend='both')
  plt.grid(True)
  ax.set_xlabel('$\mathbf{x(km)}$')
  ax.set_ylabel('$\mathbf{y(km)}$')
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  startx, endx = ax.get_xlim()
  starty, endy = ax.get_ylim()
  ax.xaxis.set_ticks(np.arange(np.floor((startx+1e-5)/(dm.Lx/1000))*dm.Lx/1000,
    np.ceil((endx-1e-5)/(dm.Lx/1000))*dm.Lx/1000+dm.Lx/1000, dm.Lx/1000))
  ax.yaxis.set_ticks(np.arange(np.floor((starty+1e-5)/(dm.Ly/1000))*dm.Ly/1000,
  np.ceil((endy-0.01)/(dm.Ly/1000))*dm.Ly/1000+dm.Ly/1000, dm.Ly/1000))

  #- colorbar
  cbar_ax = fig.add_axes([0.9, 0.3, 0.03, 0.4])
  cbar = plt.colorbar(caxc,
      ticks = np.arange(level_bd[0], level_bd[1]+1, 1), cax=cbar_ax)
  yticklabels = np.power(10.0, np.arange(level_bd[0], level_bd[1]+1, 1))
  formatting_fun = np.vectorize(lambda f: format(f, '4.0e'))
  cbar.ax.set_yticklabels(formatting_fun(yticklabels))
  cbar.ax.set_xlabel(r'$\mathbf{kg/m^3}$')
  cbar.ax.xaxis.set_label_position('top')

  #- save file
  fn = 'srfCon{0:s}_endless.png'.format(id_con[ind])
  plt.savefig(fn)
### END Main Body ###
