import numpy as np
import matplotlib.pyplot as plt

A = np.array([0, 0])
B = np.array([1, 0])
uA = 0
uB = 0

C = np.array([-0.5, 1])

P = lambda t : A*(1-t) + B*t
interp = lambda t: uA*(1-t) + uB*t
dist = lambda x : np.sqrt(np.dot(C-x, C-x))
deriv = lambda x : uB - uA + np.dot(C-P(x), A-B)/dist(P(x))

t = np.linspace(0, 1, 100)
d1 = [interp(x) for x in t]
d2 = [dist(P(x)) for x in t]
dd = [deriv(x) for x in t]


plt.figure()
plt.plot(t, d1, label="interp")
plt.plot(t, d2, label="dist")
plt.plot(t, dd)
plt.ylim([-2, 2])
plt.legend(loc=0)
plt.show()
