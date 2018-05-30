import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

vel = np.load("data.npy");
dJ = np.load("dJ.npy")

plt.matshow(np.transpose(vel), cmap="gray")
plt.colorbar()
plt.savefig("vel.png",format="png")

plt.matshow(np.transpose(dJ), cmap="gray")
plt.colorbar()
plt.savefig("dJ.png",format="png")