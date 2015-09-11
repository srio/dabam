import numpy as np

a=np.loadtxt("Mourad_Xslopes S5-06Oct-a_moy-t2.txt")

npoints = a.shape[1]

x = np.linspace(-0.5,0.5,npoints)

b = np.zeros((1+a.shape[0],npoints))

b[0,:] = x
j = 0
for i in range(5,11):
    j += 1
    b[j,:] = a[i-5,:]
for i in range(5):
    j += 1
    b[j,:] = a[i,:]


from matplotlib import pylab as plt
plt.figure(1)
for i in range(a.shape[0]):
      plt.plot(b[0,:],b[i+1,:])
plt.show()
