import numpy as np
import matplotlib.pyplot as plt
import dnest4.classic as dn4

data = dn4.my_loadtxt("data.txt")
posterior_sample = dn4.my_loadtxt("posterior_sample.txt")

for i in range(0, posterior_sample.shape[0]):
    plt.subplot(2, 1, 1)
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt="ko", alpha=0.5)
    plt.hold(True)
    plt.plot(data[:,0], posterior_sample[i, 0:data.shape[0]], "g",\
                            alpha=0.2, linewidth=2)
    plt.xlabel("Time")
    plt.ylabel("Brightness")

    plt.subplot(2, 1, 2)
    img = posterior_sample[i, data.shape[0]:].reshape((101, 201))
    plt.imshow(img, interpolation="nearest", cmap="viridis")
    plt.show()

