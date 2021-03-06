import numpy as np

a=np.loadtxt("Metz_P12128.dat",skiprows=14)

npoints = a.shape[1]
nprofiles = a.shape[0]

print("Found %d profiles, %d points each"%(nprofiles,npoints))

x = np.linspace(-0.1,0.1,npoints)
# 
# b = np.zeros((2+nprofiles,npoints))
# print("Shape a: ",a.shape," shape b: ",b.shape)
# 
# b[0,:] = x
# #default is central profile
# b[1,:] = a[nprofiles/2,:]
# 
# #recopy all profiles
# itot = 0
# for i in range(2,2+nprofiles):
#     print("i=",i)
#     b[i,:] = a[i-2,:]
#     itot += 1
# print("itot: ",itot)

# #dump file
# outFile = "Mourad.dat"
# #dd=numpy.concatenate( (x, f) ,axis=0).reshape(2,-1).transpose()
# dd=b.transpose()
# np.savetxt(outFile,dd)
# print ("File "+outFile+" written to disk:\n")


from matplotlib import pylab as plt
plt.figure(1)
for i in range(a.shape[0]):
      plt.plot(x,a[i,:])
plt.show()
