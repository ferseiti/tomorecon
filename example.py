import numpy
import matplotlib.pyplot as plt
import sys
    

def phantom( t, params ):

    a = params[0]
    n = params[1]

    if t==0:
        
        f = numpy.loadtxt('furu.dat')[0:692:8,0:692:8]

    elif t==1:

        # norma-2 (circulo)
        
        x = numpy.linspace(-a,a,n)    
        xx,yy = numpy.meshgrid(x,x)

        f = (xx**2 + yy**2 < 0.3**2).astype(numpy.double)

    elif t==2:

        # norma-1 (losango)
        
        x = numpy.linspace(-a,a,n)    
        xx,yy = numpy.meshgrid(x,x)

        f = ( numpy.abs(xx) + numpy.abs(yy) < 0.3 ).astype(numpy.double)

    else:

        # norma-p ( `quadrado`, p->infty)
        
        p = 10
        
        x = numpy.linspace(-a,a,n)    
        xx,yy = numpy.meshgrid(x,x)

        f = ( numpy.abs(xx)**p + numpy.abs(yy)**p < 0.3**p ).astype(numpy.double)
        
        
    return f


def radon( img, params ) :

    a = params[0]
    R = params[1]
    V = params[2]
    
    t  = numpy.linspace(-a,a, R)
    th = numpy.linspace( (numpy.pi * params[3] / 180.0), (numpy.pi * params[4] / 180.0), V, endpoint=False) 
    s  = numpy.linspace(-a,a,R)

    sino = numpy.zeros([R,V])   

    dx = 2.*a/(R-1)
    dy = dx
    
    for i in range(R):

        for j in range(V):
            
            for k in range(R):

                x = t[i] * numpy.cos(th[j]) - s[k] * numpy.sin(th[j])
                y = t[i] * numpy.sin(th[j]) + s[k] * numpy.cos(th[j])
        
                ix = ( numpy.ceil( (x + a)/dx) ).astype(numpy.int)
                iy = ( numpy.ceil( (y + a)/dy) ).astype(numpy.int)

                if ((ix > 0) & (ix < R) & (iy >0) & (iy < R)):         
                    sino[i][j] += img[ix][iy]
                    
    return sino


#

def backp ( sino, params ):

    a = params[0]
    R = params[1]
    V = params[2]
    
    th = numpy.linspace((numpy.pi * params[3] / 180.0), (numpy.pi * params[4] / 180.0), V, endpoint=False)
    
    x = numpy.linspace(-a,a, R)
    y = x

    dt = (2*a)/(R-1)

    b = numpy.zeros([R,R])
    
    for i in range(R):
        for j in range(R):

            cumsum = 0

            for k in range(V):

                t = x[i] * numpy.cos(th[k]) + y[j] * numpy.sin(th[k])

                idx = numpy.ceil((t + a)/dt).astype(numpy.int) 

                if ((idx > 0) & (idx < R)):
                    cumsum += sino[idx][k]

            b[i][j] = cumsum

    return b

#

def filter( sino, params ):

    a = params[0]
    R = params[1]
    V = params[2]
    
    th = numpy.linspace((numpy.pi * params[3] / 180.0), (numpy.pi * params[4] / 180.0), V, endpoint=False)
    t  = numpy.linspace(-a, a, R)
    
    dt = (2*a)/(R-1)
    
    wc = 1.0/(2*dt)
    
    w = numpy.linspace(-wc, wc, R)
    
    h = numpy.abs(w)
    
    G = numpy.fft.fftshift(numpy.transpose(numpy.kron(numpy.ones((V, 1)), h)))
    
    B = numpy.fft.fft(sino, axis=0)
    
    C = B * G
    
    #---
    
    D = numpy.fft.ifft(C, axis=0)
    
    return D.real 


#
def volRender(vol):
    
    from mayavi import mlab
    
    mlab.figure(bgcolor=(1, 1, 1))
    
    mlab.pipeline.volume(mlab.pipeline.scalar_field(vol), vmin=0, vmax=vol.max().max())
    
    mlab.outline()
    
    mlab.show()


def volContour(vol):
    
    from mayavi import mlab
    
    mlab.figure(bgcolor=(1, 1, 1))
    
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vol),plane_orientation='x_axes',slice_index=10,)
    
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vol),plane_orientation='y_axes',slice_index=10,)
    
    mlab.outline()
    mlab.show()

    
#######

'''
def phantom3D( params ):

    a = params[0]
    R = params[1]

    x = numpy.linspace()

'''
############


def tarugo( params  ):

    a = params[0]
    N = params[1]

    x = numpy.linspace(-1,1,N)
    xx,yy,zz = numpy.meshgrid(x,x,x)
    
    sphere = (xx**2 + zz**2 + yy**2 <  0.4**2).astype(numpy.double)

    cld = (xx**2 + yy**2 < 0.7**2).astype(numpy.double)
    
    vol = sphere + cld

    return vol    

####

params = ( 1.0, 74, 1, 90, 90)

treisd = tarugo( params)

volContour( treisd )

plt.imshow( treisd[37][:][:]  )
plt.show()

##

def projection( volume) :

    img = numpy.zeros([volume.shape[0],volume.shape[1]])
    
    for z in range(volume.shape[0]):
        
        img += volume[:][z][:]

    return img.T
    

plt.imshow( projection( treisd ) )
plt.show()


'''
f = phantom( 0, params )

s = radon( f, params )

plt.plot(s)
plt.show()

'''


'''
h = filter(s, params)

b = backp( h, params )

plt.figure(0)
plt.imshow(f)

plt.figure(1)
plt.imshow(s)

plt.figure(2)
plt.imshow(h)

plt.figure(3)
plt.imshow(b)

plt.show()
'''




