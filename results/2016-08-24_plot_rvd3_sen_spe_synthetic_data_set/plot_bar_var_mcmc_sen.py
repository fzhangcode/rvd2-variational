import numpy as np
import matplotlib.pyplot as plt


step = 5

# sensitivity
#mcmc = (0,0,0,0,1, 0,0,0.21,1,1, 0.14,1,1,1,1, 0.86,1,1,1,1)
#vi = (0,0,0,0.57,1, 0.07,0,0.29,1,1, 0.29,1,1,1,1, 1,0.93,0.93,1,1)

#specificity
mcmc = (1,1,1,1,1, 1,1,1,1,1, 1,0.99,0.98,1,1, 0.97,0.85,0.87,1,1)
vi = (1,1,1,1,1, 1,1,1,1,1, 1,0.98,0.98,1,1, 1,0.91,0.95,1,1)

fig, ax = plt.subplots()

index = np.arange(step)
bar_width = 0.1

fs = 20
rects1 = plt.bar(index, mcmc[0:step], bar_width, alpha=0.3, color='blue', label='MCMC, depths = 39')
rects2 = plt.bar(index + bar_width, vi[0:step], bar_width, alpha=0.3, color='magenta', label='Variational, depths = 39')

rects3 = plt.bar(index + 2*bar_width + 0.02, mcmc[step:2*step], bar_width, alpha=0.5, color='blue', label='MCMC, depths = 408')
rects4 = plt.bar(index + 3*bar_width + 0.02, vi[step:2*step], bar_width, alpha=0.5, color='magenta', label='Variational, depths = 408')

rects5 = plt.bar(index + 4*bar_width + 0.04, mcmc[2*step:3*step], bar_width, alpha=0.7, color='blue', label='MCMC, depths = 4129')
rects6 = plt.bar(index + 5*bar_width + 0.04, vi[2*step:3*step], bar_width, alpha=0.7, color='magenta', label='Variational, depths = 4129')

rects7 = plt.bar(index + 6*bar_width + 0.06, mcmc[3*step:4*step], bar_width, alpha=1, color='blue', label='MCMC, depths = 39')
rects8 = plt.bar(index + 7*bar_width + 0.06, vi[3*step:4*step], bar_width, alpha=1, color='magenta', label='Variational, depths = 41449')



def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.02*height,
                '%.2f' % height,
                ha='center', va='bottom', fontsize=15)

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
autolabel(rects5)
autolabel(rects6)
autolabel(rects7)
autolabel(rects8)

plt.setp(plt.gca().get_xticklabels(), fontsize=fs)
plt.setp(plt.gca().get_yticklabels(), fontsize=fs)

plt.ylim([0, 1.1])
plt.xlabel('VAF', fontsize = fs)
#plt.ylabel('Sensitivity', fontsize = fs)
plt.ylabel('Specificity', fontsize = fs)
plt.xticks(index + 5*bar_width, ('0.1%', '0.3%', '1.0%', '10.0%', '100.0%'))
plt.legend(loc = 'lower right', fontsize = fs)


plt.show()