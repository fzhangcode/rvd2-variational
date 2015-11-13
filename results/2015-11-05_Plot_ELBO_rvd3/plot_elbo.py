import numpy as np
import matplotlib.pyplot as plt

width = 3
msize = 10

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

filename = "./ELBO_10.txt"
elbo = []
with open(filename) as f:
    for line in f:
        parts = line.split() # split line into parts
        if len(parts) > 1:   # if at least 2 parts/columns
            elbo.append(parts[2]) 
         
plt.plot(range(len(elbo)), elbo, marker='^', linestyle='-', color='b', lw = width, markersize = msize, alpha = 0.7)

a1 = 0, 1, 2, 3, 4, len(elbo)-1
b1 = elbo[0], elbo[1], elbo[2], elbo[3], elbo[4], elbo[len(elbo)-1]  
for i, j in zip(a1, b1):
    ax.annotate('%s' %j, xy=(i,j), xytext=(8, -8),textcoords='offset points', fontsize = 20)

   
plt.xlabel('Iterations', fontsize = 30)
plt.ylabel('ELBO', fontsize = 30)
plt.xlim(-1, len(elbo) + 12)
plt.ylim(-100000000, 2000000)
plt.setp(plt.gca().get_xticklabels(), fontsize=30)
plt.setp(plt.gca().get_yticklabels(), fontsize=30)
plt.tight_layout()
plt.grid()

plt.savefig('converge_ELBO.png')