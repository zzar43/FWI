import matplotlib.pyplot as plt
import numpy as np

vel = np.load("vel.npy");

plt.matshow(np.transpose(vel), cmap="gray")
plt.colorbar()
plt.savefig("vel.png",format="png")
