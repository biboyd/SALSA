import numpy as np
import matplotlib.pyplot as plt

impact_param = np.load("impact_param.npy")

plt.hist(impact_param, bins=25)

plt.title("Distribution of Impact Parameter")
plt.xlabel("Impact Parameter (kpc)")
plt.ylabel("Counts")
plt.show()
