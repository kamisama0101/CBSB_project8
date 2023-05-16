import numpy as np

#heatmap_drawing
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap


ma1_40 = np.load("dyn41-40ma.npy")
ma41_190 = np.load("41-190ma.npy")
no_dyn=np.load("1-40ma.npy")
print(ma1_40[1:40,])
print(ma41_190[41:190,])
ma=np.concatenate((ma1_40[1:40,],ma41_190[41:190,]),axis=0)
#ma[1:100,1:100]=no_dyn[1:100,1:100]
print(ma)

data=pd.DataFrame(ma.T)
print(data)
#file = open("D:\\acdemic resources\\semester\\semester the sixth\\CBSB\\ICA\\项目\\第一篇文章\ma.txt", 'w')
#np.save(str(1)+"-"+str(A1)+'ma.npy',ma)

my_colors = ['blue','green', 'yellow']
my_cmap = ListedColormap(my_colors)
bounds = [0.1,0.125,0.145,0.165]
my_norm = BoundaryNorm(bounds, ncolors=len(my_colors))
plot=sns.heatmap(data,cmap=my_cmap,norm=my_norm)
#https://blog.csdn.net/ztf312/article/details/102474190
plot.invert_yaxis()
plt.xlabel("A1",size=10)
plt.ylabel("ST",size=10)
plt.title("Expection of a1/as",size=20)
plt.show()
