#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib import rcParams
rcParams["font.size"]=14


# 1. *Obtain the Galactic Stellar Mass Function (GSMF) exhibited by simulations Ref-L100N1504 and Recal-L025N0752 at redshifts $z = 0, 1, 2.01, \text{ and } 3.02$; for galaxies with masses $M_* > 10^8 M_\odot$.*
# 
#    Repeat for galaxies with total stellar mass computed within a fixed aperture of 30 pkpc.
# 
#    Compare and discuss the results.
# 
#    GSMF: indicates the numerical density of galaxies $(n)$ of a given stellar mass as a function of said mass.
# 
#    $$\text{GSMF}(M_*) = \log \left( \frac{dn}{d \log (M_*)} \right) \ [\text{cMpc}^{-3}]$$
# 
#    **Reference**: Furlong et al. 2015
# 
# 

# ## Primero lo hacemos para todo el volumen

# In[3]:


recal25_gmsf_s28 = pd.read_csv("recal25_gmsf_s28.txt", skiprows=15)
recal25_gmsf_s19 = pd.read_csv("recal25_gmsf_s19.txt", skiprows=15)
recal25_gmsf_s15 = pd.read_csv("recal25_gmsf_s15.txt", skiprows=15)
recal25_gmsf_s12 = pd.read_csv("recal25_gmsf_s12.txt", skiprows=15)


# In[4]:


refL0100_gmsf_s28 = pd.read_csv("refL0100_gmsf_s28.txt", skiprows=15)
refL0100_gmsf_s19 = pd.read_csv("refL0100_gmsf_s19.txt", skiprows=15)
refL0100_gmsf_s15 = pd.read_csv("refL0100_gmsf_s15.txt", skiprows=15)
refL0100_gmsf_s12 = pd.read_csv("refL0100_gmsf_s12.txt", skiprows=15)


# In[5]:


def calculate_gsmf(df, size):

    df['log_mass'] = np.log10(df['StellarMass'])
    
    bin_size = 0.3
    bins = np.arange(min(df['log_mass']), max(df['log_mass']) + bin_size, bin_size)
    
    counts, bin_edges = np.histogram(df['log_mass'], bins=bins, density=False)
    bin_midpoints = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    number_density = counts / (bin_size * size**3)
    
    log_number_density = np.log10(number_density)
    
    log_number_density[log_number_density == -np.inf] = -10  
    
    return bin_midpoints, log_number_density


# In[92]:


# List of filenames
filenames_recal = ["recal25_gmsf_s28.txt", "recal25_gmsf_s19.txt", "recal25_gmsf_s15.txt", "recal25_gmsf_s12.txt"]
size_recal = 25  

filenames_ref = ["refL0100_gmsf_s28.txt", "refL0100_gmsf_s19.txt", "refL0100_gmsf_s15.txt", "refL0100_gmsf_s12.txt"]
size_ref = 100

# Lists to store results
bin_midpoints_list_recal = []
log_number_density_list_recal = []

bin_midpoints_list_ref = []
log_number_density_list_ref = []


for filename in filenames_recal:
    df = pd.read_csv(filename, skiprows=15)    
    bin_midpoints, log_number_density = calculate_gsmf(df, size_recal)
    bin_midpoints_list_recal.append(bin_midpoints)
    log_number_density_list_recal.append(log_number_density)
    
for filename in filenames_ref:
    df = pd.read_csv(filename, skiprows=15)    
    bin_midpoints, log_number_density = calculate_gsmf(df, size_ref)
    bin_midpoints_list_ref.append(bin_midpoints)
    log_number_density_list_ref.append(log_number_density)


fig, ax = plt.subplots(1, 1, figsize=(5,4), dpi=200)



colors = ["#FF9C11", "#FF4365","#00D9C0","#292F36"]

# Plot each dataset with specific colors
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_recal, log_number_density_list_recal, labels_recal)):
    plt.plot(bin_midpoints, log_number_density, label=label, color=colors[i],linestyle='--')

# Plot each dataset with specific colors
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_ref, log_number_density_list_ref, labels_ref)):
    plt.plot(bin_midpoints, log_number_density, label=label,color=colors[i])
    
    

# Customize the plot to match the example
plt.xlabel(r'Log$_{10}$ M$_*$[M$_\odot$]', fontsize=12)
plt.ylabel(r'Log$_{10}$dn/dLog$_{10}$(M$_*$)[$\mathrm{cMpc}^{-3}]$', fontsize=12)
#plt.title('Galactic Stellar Mass Function', fontsize=14)



custom_legend = [
    Line2D([0], [0], color='black', lw=2),
    Line2D([0], [0], linestyle='--', color='black', lw=2),
    Line2D([0], [0], color=colors[3]),
    Line2D([0], [0], color=colors[2]),
    Line2D([0], [0], color=colors[1]),
    Line2D([0], [0], color=colors[0])
]

# Add custom legend
plt.legend(custom_legend, ['Recal-L025N0752', 'Ref-L100N1504 ', 'z=0', 'z=1', 'z=2.01', 'z=3.02'], loc='best',fontsize=8, framealpha=0)

# Apply tick parameters
ax = plt.gca()  # Get the current Axes instance
ax.tick_params(bottom=True, top=True, left=True, right=True, direction='in', size=5)
ax.tick_params(which="minor", bottom=True, top=True, left=True, right=True, direction='in', size=3)
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(1/ 10))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1 / 5))




plt.tight_layout()

# Show the plot
plt.show()


# ## Ahora para una apertura de 30kpc

# In[6]:


def calculate_gsmf_ap30(df, size):

    df['log_mass'] = np.log10(df['Mass'])
    
    bin_size = 0.3
    bins = np.arange(min(df['log_mass']), max(df['log_mass']) + bin_size, bin_size)
    
    # Calculate the histogram with density=True
    counts, bin_edges = np.histogram(df['log_mass'], bins=bins, density=False)
    bin_midpoints = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    number_density = counts / (bin_size * size**3)
    
    log_number_density = np.log10(number_density)
    
    log_number_density[log_number_density == -np.inf] = np.nan  
    
    return bin_midpoints, log_number_density


# In[7]:


recal25_ap30_gmsf_s28 = pd.read_csv("recal25_ap30_gmsf_s28.txt", skiprows=22)
recal25_ap30_gmsf_s19 = pd.read_csv("recal25_ap30_gmsf_s19.txt", skiprows=22)
recal25_ap30_gmsf_s15 = pd.read_csv("recal25_ap30_gmsf_s15.txt", skiprows=22)
recal25_ap30_gmsf_s12 = pd.read_csv("recal25_ap30_gmsf_s12.txt", skiprows=22)


# In[8]:


refL0100_ap30_gmsf_s28 = pd.read_csv("ref100_ap30_gmsf_s28.txt", skiprows=22)
refL0100_ap30_gmsf_s19 = pd.read_csv("ref100_ap30_gmsf_s19.txt", skiprows=22)
refL0100_ap30_gmsf_s15 = pd.read_csv("ref100_ap30_gmsf_s15.txt", skiprows=22)
refL0100_ap30_gmsf_s12 = pd.read_csv("ref100_ap30_gmsf_s12.txt", skiprows=22)


# In[89]:


filenames_recal = ["recal25_ap30_gmsf_s28.txt", "recal25_ap30_gmsf_s19.txt", "recal25_ap30_gmsf_s15.txt", "recal25_ap30_gmsf_s12.txt"]
size_recal = 25  

filenames_ref = ["ref100_ap30_gmsf_s28.txt", "ref100_ap30_gmsf_s19.txt", "ref100_ap30_gmsf_s15.txt", "ref100_ap30_gmsf_s12.txt"]
size_ref = 100

# Lists to store results
bin_midpoints_list_recal = []
log_number_density_list_recal = []

bin_midpoints_list_ref = []
log_number_density_list_ref = []


for filename in filenames_recal:
    df = pd.read_csv(filename, skiprows=22)    
    bin_midpoints, log_number_density = calculate_gsmf_ap30(df, size_recal)
    bin_midpoints_list_recal.append(bin_midpoints)
    log_number_density_list_recal.append(log_number_density)
    
for filename in filenames_ref:
    df = pd.read_csv(filename, skiprows=22)    
    bin_midpoints, log_number_density = calculate_gsmf_ap30(df, size_ref)
    bin_midpoints_list_ref.append(bin_midpoints)
    log_number_density_list_ref.append(log_number_density)

plt.figure(figsize=(5, 4),dpi=200)

# Labels for each dataset
labels_recal = ['Recal25 GMSF S28', 'Recal25 GMSF S19', 'Recal25 GMSF S15', 'Recal25 GMSF S12']
labels_ref = ['Ref100 GMSF S28', 'Ref100 GMSF S19', 'Ref100 GMSF S15', 'RRef100 GMSF S12']


colors = ["#FF9C11", "#FF4365","#00D9C0","#292F36"]
# Plot each dataset with specific colors
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_recal, log_number_density_list_recal, labels_recal)):
    plt.plot(bin_midpoints, log_number_density, label=label, color=colors[i],linestyle='--')

# Plot each dataset with specific colors
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_ref, log_number_density_list_ref, labels_ref)):
    plt.plot(bin_midpoints, log_number_density, label=label,color=colors[i])
    
    

# Customize the plot to match the example
plt.xlabel(r'Log$_{10}$ M$_*$[M$_\odot$]', fontsize=10)
plt.ylabel(r'Log$_{10}$dn/dLog$_{10}$(M$_*$)[$\mathrm{cMpc}^{-3}]$', fontsize=10)
plt.title('Galactic Stellar Mass Function (30 pkpc Aperture)', fontsize=14)



custom_legend = [
    Line2D([0], [0], color='black', lw=2),
    Line2D([0], [0], linestyle='--', color='black', lw=2),
    Line2D([0], [0], color=colors[0]),
    Line2D([0], [0], color=colors[1]),
    Line2D([0], [0], color=colors[2]),
    Line2D([0], [0], color=colors[3])
]

# Add custom legend
plt.legend(custom_legend, ['Recal-L025N0752', 'Ref-L100N1504 ', 'z=0', 'z=1', 'z=2.01', 'z=3.02'], loc='best',fontsize=8,framealpha = 0)

# Apply tick parameters
ax = plt.gca()  # Get the current Axes instance
ax.tick_params(bottom=True, top=True, left=True, right=True, direction='in', size=5)
ax.tick_params(which="minor", bottom=True, top=True, left=True, right=True, direction='in', size=3)
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(1 / 5))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1 / 5))



plt.tight_layout()

# Show the plot
plt.show()


# 
# # Probando unirlos

# In[9]:


# List of filenames
filenames_recal = ["recal25_gmsf_s28.txt", "recal25_gmsf_s19.txt", "recal25_gmsf_s15.txt", "recal25_gmsf_s12.txt"]
size_recal = 25

filenames_ref = ["refL0100_gmsf_s28.txt", "refL0100_gmsf_s19.txt", "refL0100_gmsf_s15.txt", "refL0100_gmsf_s12.txt"]
size_ref = 100

filenames_recal_30 = ["recal25_ap30_gmsf_s28.txt", "recal25_ap30_gmsf_s19.txt", "recal25_ap30_gmsf_s15.txt", "recal25_ap30_gmsf_s12.txt"]
size_recal_30 = 25

filenames_ref_30 = ["ref100_ap30_gmsf_s28.txt", "ref100_ap30_gmsf_s19.txt", "ref100_ap30_gmsf_s15.txt", "ref100_ap30_gmsf_s12.txt"]
size_ref_30 = 100

# Lists to store results
bin_midpoints_list_recal = []
log_number_density_list_recal = []

bin_midpoints_list_ref = []
log_number_density_list_ref = []

bin_midpoints_list_recal_30 = []
log_number_density_list_recal_30 = []

bin_midpoints_list_ref_30 = []
log_number_density_list_ref_30 = []

labels_recal = ['Recal z=0', 'Recal z=1', 'Recal z=2.01', 'Recal z=3.02']
labels_ref = ['Ref z=0', 'Ref z=1', 'Ref z=2.01', 'Ref z=3.02']

# Read data for recal
for filename in filenames_recal:
    df = pd.read_csv(filename, skiprows=15)
    bin_midpoints, log_number_density = calculate_gsmf(df, size_recal)
    bin_midpoints_list_recal.append(bin_midpoints)
    log_number_density_list_recal.append(log_number_density)

# Read data for ref
for filename in filenames_ref:
    df = pd.read_csv(filename, skiprows=15)
    bin_midpoints, log_number_density = calculate_gsmf(df, size_ref)
    bin_midpoints_list_ref.append(bin_midpoints)
    log_number_density_list_ref.append(log_number_density)

# Read data for recal_30
for filename in filenames_recal_30:
    df = pd.read_csv(filename, skiprows=22)
    bin_midpoints, log_number_density = calculate_gsmf_ap30(df, size_recal_30)
    bin_midpoints_list_recal_30.append(bin_midpoints)
    log_number_density_list_recal_30.append(log_number_density)

# Read data for ref_30
for filename in filenames_ref_30:
    df = pd.read_csv(filename, skiprows=22)
    bin_midpoints, log_number_density = calculate_gsmf_ap30(df, size_ref_30)
    bin_midpoints_list_ref_30.append(bin_midpoints)
    log_number_density_list_ref_30.append(log_number_density)


# In[25]:


# Create the subplots
fig, axs = plt.subplots(2, 1, figsize=(5, 8), dpi=200,sharex=True)

colors = ["#292F36", "#00D9C0", "#FF4365", "#FF9C11"]

# Plot recal datasets
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_recal, log_number_density_list_recal, labels_recal)):
    axs[0].plot(bin_midpoints, log_number_density, label=label, color=colors[i], linestyle='--')

# Plot ref datasets
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_ref, log_number_density_list_ref, labels_ref)):
    axs[0].plot(bin_midpoints, log_number_density, label=label, color=colors[i])
    
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_recal_30, log_number_density_list_recal_30, labels_recal)):
    axs[1].plot(bin_midpoints, log_number_density, label=label, color=colors[i], linestyle='--')

# Plot ref datasets
for i, (bin_midpoints, log_number_density, label) in enumerate(zip(bin_midpoints_list_ref_30, log_number_density_list_ref_30, labels_ref)):
    axs[1].plot(bin_midpoints, log_number_density, label=label, color=colors[i])
    
axs[1].set_xlabel(r'Log$_{10}$ M$_*$[M$_\odot$]', fontsize=15)
#axs[1].set_ylabel(r'Log$_{10}$dn/dLog$_{10}$(M$_*$)[$\mathrm{cMpc}^{-3}]$', fontsize=15)

axs[1].set_xlim(8,12.5)
#axs[0].sharex(axs[1])
    
for ax in axs:

    ax.set_ylabel(r'Log$_{10}$dn/dLog$_{10}$(M$_*$)[$\mathrm{cMpc}^{-3}]$', fontsize=15)
    ax.tick_params(bottom=True, top=True, left=True, right=True, direction='in', size=5,labelsize=12)
    ax.tick_params(which="minor", bottom=True, top=True, left=True, right=True, direction='in', size=3)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(1 / 3))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(1 / 3))
    
    

# Set y-axis label only for the first plot

# Add custom legend
custom_legend = [
    Line2D([0], [0], linestyle='--', color='black', lw=2),
    Line2D([0], [0], color='black', lw=2),
    Line2D([0], [0], color=colors[0]),
    Line2D([0], [0], color=colors[1]),
    Line2D([0], [0], color=colors[2]),
    Line2D([0], [0], color=colors[3])
]

fig.legend(custom_legend, ['Recal-L025N0752', 'Ref-L100N1504', 'z=0', 'z=1', 'z=2.01', 'z=3.02'], loc=(0.2, 0.1), fontsize=10, framealpha=0)

plt.tight_layout()

# Show the plot
plt.show()


# ## 3. For the simulations Ref-L0050N0752 and NoAGNL0050N0752, at redshifts $z = 0, 1, 2.01$, and $3.02$, evaluate the specific star formation rate (sSFR) as a function of the stellar mass for galaxies with $M^* > 10^8 M_\odot$ within an aperture of $30$ pkpc and sSFR $> 10^{-11}$ yr$^{-1}$ (which we will call star-forming galaxies).
# 

# In[2]:


NoAGN_ap30_fmzr_s28 = pd.read_csv("NoAGN_ap30_fmzr_s28.txt", skiprows=26)
NoAGN_ap30_fmzr_s19 = pd.read_csv("NoAGN_ap30_fmzr_s19.txt", skiprows=26)
NoAGN_ap30_fmzr_s15 = pd.read_csv("NoAGN_ap30_fmzr_s15.txt", skiprows=26)
NoAGN_ap30_fmzr_s12 = pd.read_csv("NoAGN_ap30_fmzr_s12.txt", skiprows=26)


# In[3]:


RefL50_ap30_fmzr_s28 = pd.read_csv("RefL50_ap30_fmzr_s28.txt", skiprows=29)
RefL50_ap30_fmzr_s19 = pd.read_csv("RefL50_ap30_fmzr_s19.txt", skiprows=29)
RefL50_ap30_fmzr_s15 = pd.read_csv("RefL50_ap30_fmzr_s15.txt", skiprows=29)
RefL50_ap30_fmzr_s12 = pd.read_csv("RefL50_ap30_fmzr_s12.txt", skiprows=29)


# In[11]:


def bin_sSFR_vs_StellarMass(df, num_bins=10):
    df['sSFR'] = df['SFR'] / df['Mass']
    
    log_min_mass = np.log10(df['Mass'].min())
    log_max_mass = np.log10(df['Mass'].max())
    log_bins = np.linspace(log_min_mass, log_max_mass, num_bins+1)
    
    df['log_mass_bin'] = pd.cut(np.log10(df['Mass']), bins=log_bins, labels=False)
    
    median_sSFR = df.groupby('log_mass_bin')['sSFR'].median()
    
    bin_centers = 10 ** ((log_bins[:-1] + log_bins[1:]) / 2)
    
    return bin_centers, median_sSFR

def calculate_percentiles(df, num_bins=10):
    df['sSFR'] = df['SFR'] / df['Mass']
    
    log_min_mass = np.log10(df['Mass'].min())
    log_max_mass = np.log10(df['Mass'].max())
    log_bins = np.linspace(log_min_mass, log_max_mass, num_bins+1)
    
    df['log_mass_bin'] = pd.cut(np.log10(df['Mass']), bins=log_bins, labels=False)
    
    p25_sSFR = df.groupby('log_mass_bin')['sSFR'].quantile(0.25)
    p75_sSFR = df.groupby('log_mass_bin')['sSFR'].quantile(0.75)
    
    return p25_sSFR, p75_sSFR


# In[18]:


datasets = [
    'NoAGN_ap30_fmzr_s12',
    'NoAGN_ap30_fmzr_s15',
    'NoAGN_ap30_fmzr_s19',
    'NoAGN_ap30_fmzr_s28',
    'RefL50_ap30_fmzr_s12',
    'RefL50_ap30_fmzr_s15',
    'RefL50_ap30_fmzr_s19',
    'RefL50_ap30_fmzr_s28'
]


# In[156]:


fig, ax = plt.subplots(1, 1, figsize=(5,4), dpi=200)

for dataset_name in datasets:
    df = globals()[dataset_name]  
    bin_centers, mean_sSFR = bin_sSFR_vs_StellarMass(df)
    p25_sSFR, p75_sSFR = calculate_percentiles(df)
    
    # Set linestyle based on dataset name
    linestyle = '-' if 'NoAGN' in dataset_name else '--'
    
    # Set color based on dataset name ending
    colors = {
        '12': "#292F36",
        '15': "#00D9C0",
        '19': "#FF4365",
        '28': "#FF9C11"
    }
    color = colors[dataset_name[-2:]]
    
    y_err_lower = mean_sSFR - p25_sSFR
    y_err_upper = p75_sSFR- mean_sSFR
    
    plt.plot(np.log10(bin_centers), np.log10(mean_sSFR), marker='', label=dataset_name, linestyle=linestyle, color=color)
    #plt.errorbar(bin_centers, mean_sSFR, yerr=[y_err_lower, y_err_upper], fmt='none')

plt.xlabel('$Log_{10} M_{*} \, [M_\odot]$', fontsize=14)
plt.ylabel('Median sSFR (yr$^{-1}$)', fontsize=14)




custom_legend = [
    Line2D([0], [0], color='black', lw=2),
    Line2D([0], [0], linestyle='--', color='black', lw=2),
    Line2D([0], [0], color=colores[0]),
    Line2D([0], [0], color=colores[1]),
    Line2D([0], [0], color=colores[2]),
    Line2D([0], [0], color=colores[3])
]

plt.legend(custom_legend, [ 'NoAGNL0050N0752','Ref-L0050N0752','z=0', 'z=1', 'z=2.01', 'z=3.02'], loc='best',fontsize=8,framealpha=0)

ax = plt.gca()  # Get the current Axes instance
ax.tick_params(bottom=True, top=True, left=True, right=True, direction='in', size=5)
ax.tick_params(which="minor", bottom=True, top=True, left=True, right=True, direction='in', size=3)
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(1/ 10))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1 / 5))



plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:
