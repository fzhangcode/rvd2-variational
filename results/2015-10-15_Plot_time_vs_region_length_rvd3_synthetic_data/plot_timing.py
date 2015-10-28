import numpy as np
import matplotlib.pyplot as plt

def f(t):
    return np.exp(-t) * np.cos(2*np.pi*t)

positions = [100,200,300,400]

# Time for differnt read depths
'''# 20100916_c2_p1.14_TCT.dc
time_10 = [723/60, 1194/60, 1910/60, 2557/60] #40000X
time_100 = [1265/60, 2127/60, 2933/60, 3726/60] #4000X
time_1000 = [1460/60, 2857/60, 4138/60, 5780/60] #400X
time_10000 = [359/60, 620/60, 943/60, 1294/60] #40X
'''

'''# 20100916_c3_p1.07_CGT.dc
time_10 = [160.57/60, 721.49/60, 1092.70/60, 1537.05/60]
time_100 = [141.54/60, 333.88/60, 539.88/60, 909.90/60]
time_1000 = [134.63/60, 245.12/60, 520.53/60, 675.83/60]
time_10000 = [217.16/60, 262.66/60, 541.06/60, 745.36/60 ]
'''

# 6 replicates of VCF=100.0%
time_10 = [1978.6/60, 2441.4/60, 4422.0/60, 6658.0/60]
time_100 = [1677.0/60, 2234.3/60, 3175.2/60, 4407.3/60]
time_1000 = [309.6/60, 587.8/60, 859.2/60, 1148.1/60 ]
time_10000 = [217.0/60, 375.87/60, 583.8/60, 717.3/60]

time_mcmc =[443/60, 861.7/60, 1018.2/60, 1411.1/60]

width = 3
msize = 6

plt.title('Timing for approximating variational posterior distribution', fontsize = 15)
plt.plot(positions, time_10000, marker='o', linestyle='-', color='g', lw = width, markersize = msize, label='27X Var') #'53X')
plt.plot(positions, time_1000, marker='o', linestyle='-', color='b', lw = width, markersize = msize, label='298X Var')  #'535X')
plt.plot(positions, time_100, marker='o', linestyle='-', color='r', lw = width, markersize = msize, label='3089X Var')  #'5584X')
plt.plot(positions, time_10, marker='o', linestyle='-', color='c', lw = width, markersize = msize, label='30590X Var')  #'55489X')  

plt.plot(positions, time_mcmc, marker='o', linestyle='--', color='k', lw = width, markersize = msize, label='MCMC')  #'55489X')  

plt.legend(loc='best')
plt.xlabel('Length of region of interest (positions)', fontsize = 15)
plt.ylabel('Time for converging (minutes)', fontsize = 15)
plt.xlim(80, 420)
#yticks = np.arange(0, 101, 10)
#plt.yticks(yticks)

#plt.show()
#plt.yscale('log')
plt.savefig('timing_var_mcmc.png')