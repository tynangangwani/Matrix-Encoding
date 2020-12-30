import numpy as np

mean=175
std=10
a=np.random.normal(mean, std, 10)
b=[round(i) for i in a]
print(b)
