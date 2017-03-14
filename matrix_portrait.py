import matplotlib.pyplot as plt
import numpy as np

a = np.genfromtxt('A.txt')
plt.imshow(a, cmap='hot', interpolation='nearest')
plt.show()
