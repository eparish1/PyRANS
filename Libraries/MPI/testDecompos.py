from pylab import *
nmx = 2
nmy = 2
size = 4
overlap = 1
x = linspace(0,1,113)
y = linspace(0,1,127)
z = linspace(0,1,1)
y,x,z = meshgrid(y,x,z)
nx,ny,nz = shape(x)
grid_indexG = zeros((size,3,2,2))
grid_index = zeros((3,2))
xrange = zeros(nmx)
xrangeF = zeros(nmx)
xrange[0] = int(nx/nmx) + overlap
xrangeF[0] = int(nx/nmx) 
for i in range(1,nmx-1):
  xrange[i] = int(nx/nmx) + overlap*2
  xrangeF[i] = int(nx/nmx)
xrange[-1] = nx - (sum(xrange[0:-1]) - (nmx-2)*overlap*2 -overlap) + overlap
xrangeF[-1] = nx - (sum(xrangeF[0:-1])) 
print(sum(xrange[:]) - (nmx-2)*overlap*2 -overlap*2)
print(sum(xrangeF))
yrange = zeros(nmy)
yrangeF = zeros(nmy)
yrange[0] = int(ny/nmy) + overlap
yrangeF[0] = int(ny/nmy) 
for i in range(1,nmy-1):
    yrange[i] = int(ny/nmy) + overlap*2
    yrangeF[i] = int(ny/nmy) 
yrange[-1] = ny - (sum(yrange[0:-1]) - (nmy-2)*overlap*2 -overlap) + overlap
yrangeF[-1] = ny - ( sum(yrangeF[0:-1]) )
print(sum(yrange[:]) - (nmy-2)*overlap*2 -overlap*2)
print(sum(yrangeF[:]) )

xranges = zeros(nmx)
xranges[:] = xrange[:]
for j in range(1,nmy):
    xrange = append(xrange,xranges[:])


grid_indexG[0,0,0,0] = 0
grid_indexG[0,1,0,0] = 0
grid_indexG[0,0,1,0] = xrange[0]
grid_indexG[0,1,1,0] = yrange[0]
y_ind = 0.
y_ind2 = 0.
for j in range(0,nmy):
    x_ind = 0.
    x_ind2 = 0.
    for i in range(0,nmx):
      grid_indexG[j*nmx + i,0,0,0] = x_ind
      grid_indexG[j*nmx + i,0,1,0] = x_ind + xrange[i]
      grid_indexG[j*nmx + i,1,0,0] = y_ind
      grid_indexG[j*nmx + i,1,1,0] = y_ind + yrange[j]
      grid_indexG[j*nmx + i,0,0,1] = x_ind2
      grid_indexG[j*nmx + i,0,1,1] = x_ind2 + xrangeF[i]
      grid_indexG[j*nmx + i,1,0,1] = y_ind2
      grid_indexG[j*nmx + i,1,1,1] = y_ind2 + yrangeF[j]
      x_ind += xrange[i] - 2.*overlap
      x_ind2 += xrangeF[i]
    y_ind += yrange[j] - 2.*overlap  
    y_ind2 += yrangeF[j]

print(grid_indexG[:,1,0,0])
print(grid_indexG[:,1,1,0])
#print(xrangeF)
#for j in range(0,nmy):
    #for i in range(0,nmx):
#   grid_indexG[j*nmx+i,0,0] = 0. + i*xrange[i]
#        grid_indexG[j*nmx+i,0,1] = grid_indexG[j*nmx+i,0,0] + xrange[i]
#        grid_indexG[j*nmx+i,1,0] = 0. + j*yrange[j]
#        grid_indexG[j*nmx+i,1,1] = grid_indexG[j*nmx+i,1,0] + yrange[j]
