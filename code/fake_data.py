import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

# Font tweaks
plt.rc("font", size=16, family="serif", serif="Computer Sans")
plt.rc("text", usetex=True)

def evaluate_star(x, y, limb_darkening_coefficient=1.0):
    """
    Evaluate surface brightness profile
    """
    rsq = x**2 + y**2
    which = rsq < 1.0

    star = np.zeros(x.shape)
    star[which] = 1.0 - limb_darkening_coefficient\
                            *(1.0 - np.sqrt(1.0 - rsq[which]))
    star /= star.sum()

    return star


# Set RNG seed
rng.seed(0)

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
data[:,1] = 15.0 - 2.5*np.log10(Y) + data[:,2]*rng.randn(len(t))

np.savetxt("easy_data.txt", data)

plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt="ko")
plt.gca().invert_yaxis()
plt.xlabel("Time")
plt.ylabel("Magnitude")
plt.title("`Easy' simulated data")
plt.show()

# Create harder data
def face(xx, yy):
    # The two eyes
    rsq = (xx + 0.3)**2 + (yy - 0.5)**2
    result = np.exp(-0.5*rsq/0.05**2)
    rsq = (xx - 0.3)**2 + (yy - 0.5)**2
    result += np.exp(-0.5*rsq/0.05**2)

    # The mouth
    result += np.exp(-(yy - xx**2)**2/0.05**2)*np.exp(-0.5*xx**2/0.2**2)
    return result

plt.imshow(np.exp(-face(x, y)), cmap="viridis", interpolation="nearest")
plt.show()

t = np.linspace(0.0, 10.0, 101)
Y = np.empty(len(t))
for i in range(0, len(t)):
    # Face position
    x_face = x0 + v*t[i]
    img = evaluate_star(x, y)*np.exp(-face(x - x_face, y))
    Y[i] = img.sum()
    print("{k}/{n}".format(k=i+1, n=len(t)))

data = np.empty((len(t), 3))
data[:,0] = t
data[:,2] = 0.001
data[:,1] = 15.0 - 2.5*np.log10(Y) + data[:,2]*rng.randn(len(t))

np.savetxt("hard_data.txt", data)

plt.clf()
plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], fmt="ko")
plt.gca().invert_yaxis()
plt.xlabel("Time")
plt.ylabel("Magnitude")
plt.title("`Hard' simulated data")
plt.show()


# Compute "Jeffreys Prior" for limb darkening parameter
ld = np.linspace(0.01, 0.99, 100)
h = 0.001
divergences = []
for _ld in ld:
    star0 = evaluate_star(x, y, _ld) + 1E-300
    star1 = evaluate_star(x, y, _ld + h) + 1E-300

    divergence = np.sum(star0 - star1 + star1*np.log(star1/star0))
    divergences.append(divergence)
divergences = np.array(divergences)
jeffreys_prior = np.sqrt(divergences)
jeffreys_prior /= jeffreys_prior.sum()

def badness(params):
    """
    Merit function for constructing an analytic approximation to the
    Jeffreys prior.
    """
    a, b, c = params
    approx = a + b*ld**c
    approx /= approx.sum()
    if np.any(approx < 0.0):
        return 1E100
    return np.sum(jeffreys_prior*np.log(jeffreys_prior/approx + 1E-300))

import scipy.optimize
result = scipy.optimize.minimize(badness, np.array([1.0, 1.0, 1.0]))
a, b, c = result["x"]
approx = a + b*ld**c
approx /= approx.sum()
print(a, b, c, badness([a, b, c]))

plt.plot(ld, jeffreys_prior, "ko-", ld, approx, "g")
plt.xlabel("Limb darkening parameter")
plt.ylabel("Jeffreys Prior")
plt.ylim(0.0)
plt.show()

