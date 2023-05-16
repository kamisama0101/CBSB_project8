import numpy as np
from math import factorial
from functools import reduce
from decimal import Decimal 
from decimal import getcontext 
import time

#for heatmap
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def pro_a1_true_value(A1, A2, K1, K2, ST):
    # initialize a empty matrix to store all the probabilities (comp1_comp2_comp3) of states and a1/as
    comp1_comp2_comp3 = np.full((A1+1, A2+1), Decimal(0))
    comp1_comp2_comp3[0,0]=1
    a1_in_as = np.full((A1+1, A2+1), Decimal(0))
    # calculate the components under different conditions
    for a1 in range(A1 + 1):
        for a2 in range(A2 + 1):
            as_ = a1 + a2 #bound targets
            if as_ < ST:
                if as_ >= 1:
                    a1_in_as[a1,a2]=Decimal(a1/as_)
                    # define the previous comp1_comp2_comp3[?，?] to calculate current comp1_comp2_comp3[a1,a2]
                    if a2==0:
                        comp1_comp2_comp3[a1,a2]=comp1_comp2_comp3[a1-1,a2]*Decimal((A1-a1+1)*(ST-as_+1)*K1/a1)
                    else:    
                        comp1_comp2_comp3[a1,a2]=comp1_comp2_comp3[a1,a2-1]*Decimal((A2-a2+1)*(ST-as_+1)*K2/a2)
            else:
                break
    #comp1_comp2_comp3 = np.array(comp1) * np.array(comp2) * np.array(comp3)
    comp1_comp2_comp3 = np.array(comp1_comp2_comp3)
    Z = np.sum(comp1_comp2_comp3)
    P_a = comp1_comp2_comp3 / Z
    print("结束啦")
    #print(comp1_comp2_comp3)
    print(A1,ST)
    return np.sum(P_a * a1_in_as)

#python did more accurate computation, but I don't know the result from R is not accurate

#for figure 2 in the paper
# A1=200;ST=1000;A2= 10*A1;K1 = 0.002;K2 = 0.001
# all_time=0
# T_start=time.time()
# ma=np.zeros((A1,ST))
# for bar in range(0,A1,10):
#     for a1 in range(bar+1,bar+11):
#         a2=10*a1
#         for sT in range(0,ST):
#             start = time.time()
#             ma[a1,sT]=pro_a1_true_value(a1,a2,K1,K2,sT)
#             end = time.time()
#             all_time+=end - start
#             print("(",a1,",",sT,")",":",all_time, 's',ma[a1,sT])
#     np.save('dyn'+str(41)+"-"+str(bar+10)+'ma.npy',ma)
# T_end=time.time()
# print(T_end-T_start, 's')

# data=pd.DataFrame(ma)
# print(data)



#test some cases
start = time.time()
A1=10;ST=100;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=15;ST=100;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=15;ST=20;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=15;ST=999;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=100;ST=100;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=200;ST=999;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=100;ST=999;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=10;ST=10000;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=10;ST=1000000;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=1000;ST=1000000;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')

start = time.time()
A1=10000;ST=100;A2= 10*A1;K1 = 0.002;K2 = 0.001
print(pro_a1_true_value(A1,A2,K1,K2,ST))
end = time.time()
print(end - start, 's')
