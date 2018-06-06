import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

vel = np.load("vel.npy")
# dJ = np.load("dJ.npy")
# vel_true = np.load("vel_true.npy")
# source_coor = np.load("source_coor.npy")
# receiver_coor = np.load("receiver_coor.npy")

plt.matshow(np.transpose(vel))
plt.colorbar()
# plt.scatter(source_coor[:,0], source_coor[:,1])
# plt.scatter(receiver_coor[:,0], receiver_coor[:,1], alpha=0.2)
plt.savefig("vel.png", format="png")

# plt.matshow(np.transpose(vel), cmap="gray")
# plt.colorbar()
# plt.savefig("vel.png",format="png")

# plt.matshow(np.transpose(dJ), cmap="gray")
# plt.colorbar()
# plt.savefig("dJ.png",format="png")