
# coding: utf-8

# In[4]:

import numpy
import matplotlib.pyplot as plt
import math


# In[5]:

def dmoon(y):
    Me = 600/607.0
    Mm = 7/607.0
    y1=numpy.array([0, 0, 0, 0])*1.0
    y1[0] = y[2]
    y1[1] = y[3]
    y1[2] = (2*y[3] + y[0] - (Me*(y[0] - Mm)/(((y[0] - Mm)**2 + y[1]**2)**1.5)) - (Mm*(y[0] - Me)/(((y[0] - Me)**2 + y[1]**2)**1.5))) 
    y1[3] = (-2*y[2] + y[1] - (Me*y[1]/(((y[0] - Mm)**2 + y[1]**2)**1.5)) - (Mm*y[1]/(((y[0] - Me)**2 + y[1]**2)**1.5)))
    return y1
def dfunc(y):
    Me = 600/607.0
    Mm = 7/607.0
    ra1=((y[0] + Mm)**2 + y[1]**2)**0.5
    ra2=((y[0] - Me)**2 + y[1]**2)**0.5
    uxx=-1+(Me/ra1**3)+(Mm/ra2**3)-3*(Me*(y[0]+Mm)**2)/ra1**5-3*(Mm*(y[0]-Me)**2)/ra2**5
    uyy=-1+(Me/ra1**3)+(Mm/ra2**3)-3*(Me*(y[1])**2)/ra1**5-3*(Mm*(y[1])**2)/ra2**5
    uxy=-3*y[1]*((Me*(y[0]+Mm)/ra1**5)+(Mm*(y[0]-Me)/ra2**5))
    df = numpy.identity(4)*1.0
    df[0][0]=0
    df[0][2]=1
    df[1][1]=0
    df[1][3]=1
    df[2][0]=-uxx
    df[2][1]=-uxy
    df[2][2]=0
    df[2][3]=2
    df[3][0]=-uxy
    df[3][1]=-uyy
    df[3][2]=-2
    df[3][3]=0
    return df  


# In[18]:

mu_bar= 5.589661960806324186803433350045499936607178228433273279456
nu=math.sqrt(-0.5*(mu_bar-2-math.sqrt(9*mu_bar**2-8*mu_bar)))
tau=-0.5*(nu**2+2*mu_bar+1)/nu
ax=0.001
xini=[0.84643180576+ax , 0 , 0 ,-ax*nu*tau]
t=0 
x= xini
phi = numpy.identity(4)
h=0.0001
phid=numpy.identity(4)
for ij in range(8):
    x= xini
    t=0
    phi=numpy.identity(4)
    while(x[1]*xini[3]>=0):
        k1=numpy.array(h*dmoon(x))
        k2=numpy.array(h*dmoon(x+0.5*k1))
        k3=numpy.array(h*dmoon(x+0.5*k2))
        k4=numpy.array(h*dmoon(x+k3))
        phid=numpy.dot(dfunc(x),phi)
        phi=phi+phid*h
        x=x+(k1 + 2*k2 + 2*k3 + k4)/6.00
        t=t+h
    delta_vy=x[2]*x[3]/(phi[2][3]*x[3]+phi[1][3]*x[2])   
    xini[3]-=delta_vy
print t
xini1=numpy.array(xini)
print xini1


# In[19]:

ax=2*ax
xini=[0.84643180576+ax , 0 , 0 ,-ax*nu*tau]
t=0 
x= xini
phi = numpy.identity(4)
h=0.0001
phid=numpy.identity(4)
for ij in range(8):
    x= xini
    t=0
    phi=numpy.identity(4)
    while(x[1]*xini[3]>=0):
        k1=numpy.array(h*dmoon(x))
        k2=numpy.array(h*dmoon(x+0.5*k1))
        k3=numpy.array(h*dmoon(x+0.5*k2))
        k4=numpy.array(h*dmoon(x+k3))
        phid=numpy.dot(dfunc(x),phi)
        phi=phi+phid*h
        x=x+(k1 + 2*k2 + 2*k3 + k4)/6.00
        t=t+h
    delta_vy=x[2]*x[3]/(phi[2][3]*x[3]+phi[1][3]*x[2])   
    xini[3]-=delta_vy
xini2=numpy.array(xini)
print xini2


# In[20]:

for ty in range(10):
    xini=2*xini2-xini1
    t=0 
    x= xini
    phi = numpy.identity(4)
    h=0.0001
    phid=numpy.identity(4)
    for ij in range(8):
        x= xini
        t=0
        phi=numpy.identity(4)
        while(x[1]*xini[3]>=0):
            k1=numpy.array(h*dmoon(x))
            k2=numpy.array(h*dmoon(x+0.5*k1))
            k3=numpy.array(h*dmoon(x+0.5*k2))
            k4=numpy.array(h*dmoon(x+k3))
            phid=numpy.dot(dfunc(x),phi)
            phi=phi+phid*h
            x=x+(k1 + 2*k2 + 2*k3 + k4)/6.00
            t=t+h
        delta_vy=x[2]*x[3]/(phi[2][3]*x[3]+phi[1][3]*x[2])   
        xini[3]-=delta_vy
    print xini
    xini1=xini2
    xini2=xini

print xini2
print t


# In[22]:

xmid=numpy.array([0.,0.,0.,0.])
epsilon=10**-3
h=0.0001
n = int(t/h)
print n
y=numpy.arange(8*n).reshape((2*n,4))*1.0
y[0]=xini2
phi = numpy.identity(4)*1.0
for i in range(1,2*n):
    k1=numpy.array(h*dmoon(y[i-1]))
    k2=numpy.array(h*dmoon(y[i-1]+0.5*k1))
    k3=numpy.array(h*dmoon(y[i-1]+0.5*k2))
    k4=numpy.array(h*dmoon(y[i-1]+k3))
    y[i]=y[i-1]+(k1+2*k2+2*k3+k4)/6.0
    phid=numpy.dot(dfunc(y[i-1]),phi)     
    phi=phi+phid*h
plt.plot(y[:,0],y[:,1])
lam , eigvectors= numpy.linalg.eig(phi)
eigvectors=numpy.real(eigvectors)
lam = (lam)
print lam 
print lam[0]*lam[1]
print eigvectors


# In[17]:

for i in range(1,2*n):
    k1=numpy.array(h*dmoon(y[i-1]))
    k2=numpy.array(h*dmoon(y[i-1]+0.5*k1))
    k3=numpy.array(h*dmoon(y[i-1]+0.5*k2))
    k4=numpy.array(h*dmoon(y[i-1]+k3))
    y[i]=y[i-1]+(k1+2*k2+2*k3+k4)/6.0
    phid=numpy.dot(dfunc(y[i-1]),phi)     
    phi=phi+phid*h
    '''
    if(i%500==0):
        h_mani=0.001
        n_mani=100
        y_manifold=numpy.arange(8*n_mani).reshape(2*n_mani,4)*0.0
        dis_vector=numpy.dot(phi,eigvectors[1])
        dis_vector=dis_vector/numpy.linalg.norm(dis_vector)
        y_manifold[0]=y[i]+epsilon*dis_vector
        for j in range(1, 2*n_mani):
            
            k1_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]))
            k2_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]+0.5*k1))
            k3_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]+0.5*k2))
            k4_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]+k3))
            y_manifold[j] = y_manifold[j-1]+(k1_mani+2*k2_mani+2*k3_mani+k4_mani)/6.0
            if((y_manifold[j][0]<0 and y_manifold[j][1]<0) or y_manifold[j][0]>600/607.0):
                break
        plt.plot(y_manifold[:,0],y_manifold[:,1])
        '''
    if(i%500==100):
        h_mani=0.001
        n_mani=1500
        y_manifold=numpy.arange(8*n_mani).reshape(2*n_mani,4)*0.0
        dis_vector=numpy.dot(phi,eigvectors[0])
        dis_vector=dis_vector/numpy.linalg.norm(dis_vector)
        y_manifold[0]=y[i]+epsilon*dis_vector
        for j in range(1, 2*n_mani):
            
            k1_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]))
            k2_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]+0.5*k1))
            k3_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]+0.5*k2))
            k4_mani=numpy.array(h_mani*dmoon(y_manifold[j-1]+k3))
            y_manifold[j] = y_manifold[j-1]+(k1_mani+2*k2_mani+2*k3_mani+k4_mani)/6.0
            if((y_manifold[j][0]<0 and y_manifold[j][1]<0) or y_manifold[j][0]>600/607.0):
                break
        plt.plot(y_manifold[:,0],y_manifold[:,1])


# In[ ]:



