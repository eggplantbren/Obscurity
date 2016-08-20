import numpy as np
import matplotlib.pyplot as plt
import dnest4.classic as dn4

data = dn4.my_loadtxt("data.txt")
posterior_sample = dn4.my_loadtxt("posterior_sample.txt")

for i in range(0, posterior_sample.shape[0]):
    line = posterior_sample[i, :].copy()

    plt.subplot(3, 1, 1)
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt="ko")
    plt.hold(True)
    plt.plot(data[:,0], line[0:data.shape[0]], "g",\
                            linewidth=2)
    plt.xlabel("Time")
    plt.ylabel("Brightness")
    line = line[data.shape[0]:]

    plt.subplot(3, 1, 2)
    img = line[0:len(line)//2].reshape((51, 101))
    plt.imshow(img, interpolation="nearest", cmap="viridis")
    line = line[len(line)//2:]

    plt.subplot(3, 1, 3)
    img = line.reshape((51, 101))
    plt.imshow(img, interpolation="nearest", cmap="viridis")
    plt.show()

