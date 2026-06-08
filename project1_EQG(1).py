#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pynbody
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob
import h5py
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import astropy.units as u
import astropy.constants as c
from matplotlib.ticker import MultipleLocator


# In[ ]:


# Path to Amanda's directory with the snapshots
path_snapshots = '/home/amanda.ormazbal/gadget4/examples/CollidingGalaxiesSFR/output/'
# We consider only the files that start with 'snapshot_' and end with '.hdf5'.
snapshots_files = glob.glob(f'{path_snapshots}/snapshot_*.hdf5')
# Here's a function to extract the snapshot number:
def extract_snapshot_number(snapshots_files):
    base = snapshots_files.split('/')[-1]
    number = base.replace('snapshot_', '').replace('.hdf5', '')
    return int(number)

snapshots_files_sorted = sorted(snapshots_files, key=extract_snapshot_number)
#print('\nsorted snapshots :) :', snapshots_files_sorted)


# In[ ]:


sim_list = []

for i in range(1, 39):
    #sim = pynbody.load(snapshots_files_sorted[i])
    sim = h5py.File(snapshots_files_sorted[i])
    sim_list.append(sim)    


# In[ ]:


file_path = '/home/amanda.ormazbal/gadget4/examples/CollidingGalaxiesSFR/output/snapshot_000.hdf5'
snapshot = h5py.File(file_path , "r")
x, y, z = snapshot["PartType0"]["Coordinates"][:,0], snapshot["PartType0"]["Coordinates"][:,1], snapshot["PartType0"]["Coordinates"][:,2]
x, y, z


# In[ ]:


def add_scalebar(ax, l, units='kpc', loc='lower right', size=16, color='k'): #original size=16
    scalebar = AnchoredSizeBar(
                        ax.transData,
                        l,
                        f'{l} {units}',
                        loc,
                        pad=0.5,
                        color=color,
                        frameon=False,
                        size_vertical=l / 60,
                        fontproperties=fm.FontProperties(size=size),
                        sep=6)
    ax.add_artist(scalebar)
    ax.set(yticks=[], xticks=[])


# In[11]:


import matplotlib.pyplot as plt
import imageio

scatter_plot_files = []


# Loop through each snapshot in sim_list
for i, sim in enumerate(sim_list):
    # Create a scatter plot of the positions of stars
    fig, ax = plt.subplots(1,2, figsize = (10,4), dpi = 200)
    
    ax[0].scatter(sim["PartType0"]["Coordinates"][:,0], sim["PartType0"]["Coordinates"][:,1],alpha=0.4, s=0.001,c='orange',zorder=2)
    ax[0].scatter(sim["PartType1"]["Coordinates"][:,0], sim["PartType1"]["Coordinates"][:,1],alpha=0.4, s=0.001,c='darkmagenta',zorder=0)
    ax[0].scatter(sim["PartType2"]["Coordinates"][:,0], sim["PartType2"]["Coordinates"][:,1],alpha=0.4, s=0.001,c='navy',zorder=1)
    ax[0].scatter(sim["PartType3"]["Coordinates"][:,0], sim["PartType3"]["Coordinates"][:,1],alpha=0.4, s=0.001,c='navy',zorder=1)
    ax[0].scatter(sim["PartType4"]["Coordinates"][:,0], sim["PartType4"]["Coordinates"][:,1],alpha=0.4, s=0.001,c='navy',zorder=1)
    add_scalebar(ax[0], 40, units='kpc', loc='lower right', size=16, color='k')
    
    #plt.title('Time {}'.format(str(sim.properties['time'])))
    #ax[0].set_xlabel('X [kpc]')
    #ax[0].set_ylabel('Y [kpc]')
    ax[0].set_xlim(-100, 100)  # Adjust as needed
    ax[0].set_ylim(-100, 100)  # Adjust as needed
    ax[0].text(-90, 85, 'Stars', color='blue', fontsize=10)
    ax[0].text(-90, 75, 'Gas', color='orange', fontsize=10)
    ax[0].text(-90, 65, 'Dark Matter', color='darkmagenta', fontsize=10)

    ax[1].scatter(sim["PartType0"]["Coordinates"][:,0], sim["PartType0"]["Coordinates"][:,2],alpha=0.4, s=0.001,c='orange',zorder=2)
    ax[1].scatter(sim["PartType1"]["Coordinates"][:,0], sim["PartType1"]["Coordinates"][:,2],alpha=0.4, s=0.001,c='darkmagenta',zorder=0)
    ax[1].scatter(sim["PartType2"]["Coordinates"][:,0], sim["PartType2"]["Coordinates"][:,2],alpha=0.4, s=0.001,c='navy',zorder=1)
    ax[1].scatter(sim["PartType3"]["Coordinates"][:,0], sim["PartType3"]["Coordinates"][:,2],alpha=0.4, s=0.001,c='navy',zorder=1)
    ax[1].scatter(sim["PartType4"]["Coordinates"][:,0], sim["PartType4"]["Coordinates"][:,2],alpha=0.4, s=0.001,c='navy',zorder=1)
    add_scalebar(ax[1], 40, units='kpc', loc='lower right', size=16, color='k')
    
   
    ax[1].set_xlim(-100, 100)  # Adjust as needed
    ax[1].set_ylim(-100, 100)  # Adjust as needed
    
    ax[1].text(0.70, 0.95, f'Snapshot {i+1}', color='black', fontsize=10, transform=ax[1].transAxes, verticalalignment='top')
    
    # Save the scatter plot as an image file
    filename = f'scatter_plot_{i+1}.png'
    scatter_plot_files.append(filename)
    plt.savefig(filename)
    plt.close()

# Create a GIF from the scatter plot images
images = []
for filename in scatter_plot_files:
    images.append(imageio.imread(filename))

# Save the GIF
output_gif = 'scatter_plots.gif'
imageio.mimsave(output_gif, images)

# Display the GIF
from IPython.display import Image
Image(open(output_gif, 'rb').read())


# In[7]:


sim_list_pnbody = []

for i in range(1, 39):
    sim = pynbody.load(snapshots_files_sorted[i])
    sim_list_pnbody.append(sim)    


# In[8]:


sim_list_pnbody[0].gas.loadable_keys()


# In[9]:


def temp ( u ) :
    mu = 0.61 # Mean molecular weight
    gamma = 5.0 / 3.0 # Adiabatic index
    kb = 1.38064852e-16 # Boltzmann ’ s constant in cgs
    mp = 1.67262192e-24 # Proton mass in cgs
    temperature = mu * ( gamma - 1) * mp * u / kb
    return temperature


# In[10]:


energia_interna_1 = sim_list_pnbody[1].gas["u"]* u.cm**2 / u.s**(2)
energia_interna_2 = sim_list_pnbody[9].gas["u"]* u.cm**2 / u.s**(2)
energia_interna_3 = sim_list_pnbody[12].gas["u"]* u.cm**2 / u.s**(2)
energia_interna_4 = sim_list_pnbody[25].gas["u"]* u.cm**2 / u.s**(2)
energia_interna_5 = sim_list_pnbody[28].gas["u"]* u.cm**2 / u.s**(2)
energia_interna_6 = sim_list_pnbody[-1].gas["u"]* u.cm**2 / u.s**(2)


# In[37]:


plt.axhline(-6.23,c='salmon')
plt.axvline(-3.23,c='salmon')

plt.scatter(np.log10(sim_list_pnbody[-1].gas["rho"]),np.log10(temp(energia_interna_6.value)), s=5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 38',c='lightblue')
plt.scatter(np.log10(sim_list_pnbody[28].gas["rho"]),np.log10(temp(energia_interna_5.value)), s=5, edgecolor = None,alpha= 0.7, lw= 0,label = 'Snapshot 28',c='orange')
plt.scatter(np.log10(sim_list_pnbody[25].gas["rho"]),np.log10(temp(energia_interna_4.value)), s=5, edgecolor = None,alpha= 0.5, lw= 0,label = 'Snapshot 25',c='green')
plt.scatter(np.log10(sim_list_pnbody[12].gas["rho"]),np.log10(temp(energia_interna_3.value)), s=5, edgecolor = None,alpha= 0.3, lw= 0,label = 'Snapshot 12',c='navy')
plt.scatter(np.log10(sim_list_pnbody[9].gas["rho"]),np.log10(temp(energia_interna_2.value)),  s=5, edgecolor = None,alpha= 0.2, lw= 0,label = 'Snapshot 9', c='#ffef33')
plt.scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 1', c='salmon')

plt.legend()





# In[39]:


plt.axhline(-6.23,c='salmon')
plt.axvline(-3.23,c='salmon')

plt.scatter(np.log10(sim_list_pnbody[-1].gas["rho"]),np.log10(temp(energia_interna_6.value)), s=5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 38',c='lightblue')
plt.scatter(np.log10(sim_list_pnbody[28].gas["rho"]),np.log10(temp(energia_interna_5.value)), s=5, edgecolor = None,alpha= 0.4, lw= 0,label = 'Snapshot 28',c='orange')
#plt.scatter(np.log10(sim_list_pnbody[25].gas["rho"]),np.log10(temp(energia_interna_4.value)), s=5, edgecolor = None,alpha= 0.5, lw= 0,label = 'Snapshot 25',c='green')
#plt.scatter(np.log10(sim_list_pnbody[12].gas["rho"]),np.log10(temp(energia_interna_3.value)), s=5, edgecolor = None,alpha= 0.3, lw= 0,label = 'Snapshot 12',c='navy')
#plt.scatter(np.log10(sim_list_pnbody[9].gas["rho"]),np.log10(temp(energia_interna_2.value)),  s=5, edgecolor = None,alpha= 0.2, lw= 0,label = 'Snapshot 9', c='#ffef33')
#plt.scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=5, edgecolor = None,alpha= 0.4, lw= 0,label = 'Snapshot 1', c='salmon')

plt.legend()


# In[70]:


fig, ax = plt.subplots(2,3, figsize = (10,6), dpi = 200)

ax[0,0].scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=0.5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 1', c='salmon')
ax[0,0].set_ylabel(r'$\text{Log}_{10} \; (T/K)$')
ax[0,0].xaxis.set_major_locator(MultipleLocator(5))
ax[0,0].text(0.25, 0.9, 'Snapshot 1', fontsize=10, color='black', transform=ax[0,0].transAxes, ha='center', va='center',bbox=dict(facecolor='white', alpha=0.7))


ax[0,1].scatter(np.log10(sim_list_pnbody[9].gas["rho"]),np.log10(temp(energia_interna_2.value)),  s=0.5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 9', c='navy')
ax[0,1].scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=0.5, edgecolor = None,alpha= 0.4, lw= 0,label = 'Snapshot 1', c='salmon')
ax[0,1].set_xlim(-15,0)
ax[0,1].set_ylim(-9,-3)
ax[0,1].text(0.25, 0.9, 'Snapshot 9', fontsize=10, color='black', transform=ax[0,1].transAxes, ha='center', va='center',bbox=dict(facecolor='white', alpha=0.7))

ax[0,2].scatter(np.log10(sim_list_pnbody[12].gas["rho"]),np.log10(temp(energia_interna_3.value)),  s=0.5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 12', c='navy')
ax[0,2].scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=0.5, edgecolor = None,alpha= 0.4, lw= 0,label = 'Snapshot 1', c='salmon')
ax[0,2].set_xlim(-15,0)
ax[0,2].set_ylim(-9,-3)
ax[0,2].text(0.25, 0.9, 'Snapshot 12', fontsize=10, color='black', transform=ax[0,2].transAxes, ha='center', va='center',bbox=dict(facecolor='white', alpha=0.7))

ax[1,0].scatter(np.log10(sim_list_pnbody[25].gas["rho"]),np.log10(temp(energia_interna_4.value)),  s=0.5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 25', c='navy')
ax[1,0].scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=0.5, edgecolor = None,alpha= 0.4, lw= 0,label = 'Snapshot 1', c='salmon')
ax[1,0].set_xlim(-15,0)
ax[1,0].set_ylim(-9,-3)
ax[1,0].set_xlabel(r'$log_{10} \; (\rho)$')
ax[1,0].set_ylabel(r'$log_{10} \; (T/K)$')
ax[1,0].text(0.25, 0.9, 'Snapshot 25', fontsize=10, color='black', transform=ax[1,0].transAxes, ha='center', va='center',bbox=dict(facecolor='white', alpha=0.7))
ax[1,0].xaxis.set_major_locator(MultipleLocator(1))

ax[1,1].scatter(np.log10(sim_list_pnbody[28].gas["rho"]),np.log10(temp(energia_interna_5.value)),  s=0.5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 28', c='navy')
ax[1,1].scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=0.5, edgecolor = None,alpha= 0.4, lw= 0,label = 'Snapshot 1', c='salmon')
ax[1,1].set_xlim(-15,0)
ax[1,1].set_ylim(-9,-3)
ax[1,1].set_xlabel(r'$log_{10} \; (\rho)$')
ax[1,1].text(0.25, 0.9, 'Snapshot 28', fontsize=10, color='black', transform=ax[1,1].transAxes, ha='center', va='center',bbox=dict(facecolor='white', alpha=0.7))

ax[1,2].scatter(np.log10(sim_list_pnbody[-1].gas["rho"]),np.log10(temp(energia_interna_6.value)),  s=0.5, edgecolor = None,alpha= 1, lw= 0,label = 'Snapshot 38', c='navy')
ax[1,2].scatter(np.log10(sim_list_pnbody[1].gas["rho"]),np.log10(temp(energia_interna_1.value)),  s=0.5, edgecolor = None,alpha= 0.4, lw= 0,label = 'Snapshot 1', c='salmon')
ax[1,2].set_xlim(-15,0)
ax[1,2].set_ylim(-9,-3)
ax[1,2].set_xlabel(r'$log_{10} \; (\rho)$')
ax[1,2].text(0.25, 0.9, 'Snapshot 38', fontsize=10, color='black', transform=ax[1,2].transAxes, ha='center', va='center',bbox=dict(facecolor='white', alpha=0.7))

ax = np.ravel(ax)

for i in range(len(ax)):
    
    ax[i].xaxis.set_major_locator(MultipleLocator(5))
    ax[i].yaxis.set_major_locator(MultipleLocator(1))
    ax[i].tick_params(bottom = True, top = False, left = True, right = False, direction = 'out', size = 3.5, width = 1.8)
    ax[i].tick_params(which ="minor",bottom = False, top = False, left = False, right = False, direction = 'in', size = 3)
    for axis in ['top','bottom','left','right']:
        ax[i].spines[axis].set_linewidth(1.7)



# In[19]:


pynbody.plot.gas.rho_T(sim_list_pnbody[-1])


# In[ ]:




