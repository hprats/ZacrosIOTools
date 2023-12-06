import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 3, 5, 8, 12, 13, 16, 17, 20, 22, 24])
t = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

a = np.array([x, t])

print(np.diff(a))
print(np.gradient(x, t))
