import numpy as np
import os
import matplotlib.pyplot as plt

data_dir = '../data/'
fig_dir = '../figs/'
os.makedirs(fig_dir,exist_ok=True)


f = open(data_dir+'nd.dac','r')
nd = int(f.read().split()[0])
f.close()

f = open(data_dir+'params.dac','r')
params = f.read().split()
f.close()

nx = int(params[0])
mtype = int(params[1])
xmax = float(params[2])
xmin = float(params[3])
gm = float(params[4])

dx = (xmax - xmin)/nx
x = np.linspace(xmin + 0.5*dx, xmax - 0.5*dx,nx)

endian = '<'

plt.close('all')
plt.clf()
fig = plt.figure(num=100,figsize=(5,6))

xm = 0.5*(xmax + xmin)
dd = 0.04

n0 = 0
n1 = 1
n1 = nd+1
for n in range(n0,n1):
    print(n)
    f = open(data_dir+'t.dac.'+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,endian+'d',1)
    f.close()
    t = t.reshape(1,order='F')[0]

    gc = xm + t 
    qqc = np.exp(-((x-gc)/dd)**2)   + np.exp(-((x-gc+xmax)/dd)**2) 

    dtype = np.dtype([("qq",endian+str(nx*mtype)+"d")])
    f = open(data_dir+'qq.dac.'+'{0:08d}'.format(n),"rb")
    qq = np.fromfile(f,dtype=dtype,count=1)
    qq = qq["qq"].reshape((nx,mtype),order='F')
    ro = qq[:,0]
    vx = qq[:,1]
    vy = qq[:,2]
    vz = qq[:,3]
    bx = qq[:,4]
    by = qq[:,5]
    bz = qq[:,6]
    se = qq[:,7]

    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.plot(x,qq[:,0],'.-')
    #ax1.plot(x,np.flip(qq[:,0])-1,'.',color='red')
    rr = 1.0
    cv = rr/(gm - 1)
    pr = qq[:,0]**gm*exp(qq[:,7]/cv)
    ax2.plot(x,pr,'.-')
    #ax2.plot(x,pr,'.-')
    #ax.set_ylim(-0.01,0.01)
    plt.pause(0.01)

    if(n == n0):
        fig.tight_layout()
    fig.savefig(fig_dir+'py'+'{0:04d}'.format(n)+'.png')

    if(n != n1-1):
        plt.clf()

        
