import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

rng.seed(0)

N = 101
t = np.linspace(-1.0, 1.0, N)
y = 1.0 - 0.1*np.exp(-0.5*(t/0.1)**2)
sig = 0.01*np.ones(N)
y += sig*rng.randn(N)

data = np.vstack([t, y, sig]).T
np.savetxt("data.txt", data)

plt.plot(t, y, "ko")
plt.show()

