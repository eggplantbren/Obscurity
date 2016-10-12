import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

rng.seed(0)

def evaluate_star(x, y, limb_darkening_coefficient=1.0):
    """
    Evaluate surface brightness profile
    """
    rsq = x**2 + y**2

    limb_darkening_coefficient = 1.0
    which = rsq < 1.0

    star = np.zeros(x.shape)
    star[which] = 1.0 - limb_darkening_coefficient\
                            *(1.0 - np.sqrt(1.0 - rsq[which]))
    star /= star.sum()

    return star


# Coordinate grid
x = np.linspace(-5, 5, 1001)
y = np.linspace(-5, 5, 1001)
[x, y] = np.meshgrid(x, y)
y = y[::-1, :]

# Initial position and speed of a blob
x0, v = -3.0, 0.5

# Blob radius
a = 0.3

# Create the data
t = np.linspace(0.0, 10.0, 101)
Y = np.empty(len(t))
for i in range(0, len(t)):
    # Blob position
    x_blob = x0 + v*t[i]

    # Squared distance of pixels from blob center
    rsq_blob = (x - x_blob)**2 + y**2
    img = evaluate_star(x, y)*(rsq_blob > a**2)

    Y[i] = img.sum()
    print("{k}/{n}".format(k=i+1, n=len(t)))

data = np.empty((len(t), 3))
data[:,0] = t
data[:,2] = 0.01
data[:,1] = Y + data[:,2]*rng.randn(len(t))

np.savetxt("data.txt", data)

plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt="ko")
plt.show()



# Compute 
