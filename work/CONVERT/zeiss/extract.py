import numpy

a = numpy.loadtxt("/home/srio/OASYS1.2/dabam/work/CONVERT/zeiss/ZEISS 19-0057_VFM.asc",skiprows=1)

print(a.shape)

y = a[:,0]
x = a[:,1]
z = a[:,2]

from srxraylib.plot.gol import plot, plot_image

i0 = numpy.argwhere(x == 0)
ny = len(i0)
nx = z.size // ny

yy = numpy.array(y[i0])
yy = yy.reshape(-1)
xx = numpy.unique(x)

zz = numpy.zeros((nx,ny))

ij = -1

for i in range(nx):
    for j in range(ny):
        ij += 1
        zz[i,j] = z[ij]

print(xx.shape,yy.shape,zz.shape, xx.size*yy.size)


zcentral = z[i0]
plot(yy, zcentral,
     yy, zz[nx//2,:])
plot_image(zz, xx, yy, aspect='auto' )

from oasys.util.oasys_util import write_surface_file
write_surface_file(1e-3*zz.T,1e-3*xx,1e-3*yy,file_name="ZEISS_19-0057_VFM.h5",overwrite=True)
print("File written to disk: ZEISS_19-0057_VFM.h5")