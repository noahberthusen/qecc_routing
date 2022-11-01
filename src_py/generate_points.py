from mec import make_circle
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt

N = 20
k = 7

all_points = np.array([(x,y) for x in range(N) for y in range(N)], dtype="i,i")
points = random.sample(list(np.arange(len(all_points))), k)

def in_circle(point, circle):
    x, y = point
    cx, cy, r = circle
    return np.sqrt((x - cx)**2 + (y - cy)**2) <= r

# gens = pd.DataFrame()
full_arr = []
rs = np.linspace(1, np.sqrt(2)*N, 20)

for r in rs:
    for point in all_points:
        for i in range(5):
            arr = []
            # points = random.sample(list(np.arange(len(all_points))), k)
            cx, cy = point
            # r = random.uniform(1, np.sqrt(2)*N)
            
            in_points = all_points[[in_circle(point, (cx, cy, r)) for point in all_points]]

            if (len(in_points) > k):
                points = random.sample(list(np.arange(len(in_points))), k)
                circle = make_circle(in_points[points])

                arr.append(list(in_points[points]))
                arr += [*circle]
                full_arr.append(arr)

gens = pd.DataFrame(full_arr, columns=["points", "cx", "cy", "r"])
print(len(gens))
gens.to_csv(f"{N}_{k}_gens.csv", index=False)
print(gens.head())
plt.hist(gens['r'])
plt.show()
