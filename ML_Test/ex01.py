import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


D = np.loadtxt("dataQuadReg2D.txt")
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(D[:,0],D[:,1],D[:,2], "ro")
plt.show()
