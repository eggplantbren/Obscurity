import numpy as np
import matplotlib.pyplot as plt
import os
import dnest4.classic as dn4

data = dn4.my_loadtxt("data.txt")
posterior_sample = dn4.my_loadtxt("posterior_sample.txt")

print("WARNING! This will delete\
 movie.mkv and the Frames/ directory, if these exist.")
ch = input("Continue? y/n: ")
if ch != "y" and ch != "Y":
    exit()

os.system("rm -rf Frames/ movie.mkv")
os.mkdir("Frames")

for i in range(0, posterior_sample.shape[0]):
    line = posterior_sample[i, :].copy()
    line = line[5:]

    plt.subplot(3, 1, 1)
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt="ko")
    plt.hold(True)
    plt.plot(data[:,0], line[0:data.shape[0]], "g",\
                            linewidth=2)
    plt.gca().invert_yaxis()
    plt.xlabel("Time")
    plt.ylabel("Brightness")
    line = line[data.shape[0]:]

    plt.subplot(3, 1, 2)
    img = line[0:len(line)//2].reshape((101, 201))
    plt.imshow(img, interpolation="nearest", cmap="viridis")
    line = line[len(line)//2:]

    plt.subplot(3, 1, 3)
    img = line.reshape((101, 201))
    plt.imshow(img, interpolation="nearest", cmap="viridis")

    plt.savefig("Frames/" + "%0.6d"%(i+1) + ".png", bbox_inches="tight")
    print("Processed frame {k}/{n}.".format(k=(i+1),
                        n=posterior_sample.shape[0]))

# Make a movie with ffmpeg
os.system("ffmpeg -r 10 -i Frames/%06d.png -c:v h264 -b:v 4192k movie.mkv")

