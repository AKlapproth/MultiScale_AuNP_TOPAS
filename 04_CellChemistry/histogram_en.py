import numpy as np
import matplotlib.pyplot as plt
from topas2numpy import read_ntuple

PhSp_name = "ChemicalHistories.phsp"
x = read_ntuple(PhSp_name)
species_label = 'Particle Type Flag'
species = x[species_label]
print('start histogram')

n, bins, patches = plt.hist(species, 6)
plt.xlabel(species_label)
plt.ylabel('Histories')
# plt.title('Photon y')
# plt.savefig('Photon y')
plt.title('Chemical Species')
plt.savefig('Chemical Species')
plt.show()


# arr = plt.hist(species, 6)
