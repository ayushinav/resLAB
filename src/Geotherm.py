
# Importing libraries
import numpy as np;
import pandas as pd;
import matplotlib.pyplot as maps;

def geotherm(Qs, Ts, A, k, z):
    """
    Function for estimating geotherm for an n-layered earth
    Qs= Surface Heat Flow
    Ts= Surface Temperature
    A= vector of heat production rate
    k= vector of thermal conductivity
    """
    
    n= len(z); # number of layers
    Q= Qs;
    T= Ts;
    z_1= 0;
    depth= np.array([]);
    gtherm= np.array([]);
    for i in range(0, n):
        print(Q*1000)
        Z= np.arange(z_1, z[i]+1, 1);
        
        temp= -A[i]/2/k[i]* Z*Z+ (Q+ A[i]*z_1)/k[i]* Z + (T- Q/k[i]*z_1- A[i]/2/k[i]* z_1**2);
        gtherm= np.hstack((gtherm, temp));
        depth= np.hstack((depth, Z));
        Q= Q+ A[i]*(z_1- z[i]);
        T= temp[-1]; 
        z_1= z[i];
        print(Q, "\t", temp[-1])
    
    print("Temperature at lowermost layer is %d K" %(T));
    print("Heat flow estimated at lowermost interface is %1.5f Wm-2 K-1" %(Q));
    return [gtherm, depth];

def plot_geotherm(ax, z, data_z, T):
    """
    Function to plot geotherm
    z= vector of [heat producing crust, Moho depth, depth to LAB]
    data_z= depth vector in km
    T= Temperature in K
    """
    ax.plot(T, data_z, label= "Geotherm");
    ax.set_ylim([z[-1]/1000+ 50, 0]);
    #ax.set_xlim([Ts, Ts+3000])
    ax.grid('on');
    moho= z[1]/1000*np.ones_like(T);
    LAB= z[2]/1000*np.ones_like(T);
    range_rho= np.logspace(1, 5, len(data_z))
    ax.fill_between(T, moho, 0*LAB, color= 'green', alpha= 0.3, label= 'crust');
    ax.fill_between(T, moho, LAB, color= 'yellow', alpha= 0.25, label= 'Lithospheric mantle');
    ax.fill_between(T, LAB, 300*LAB/LAB, color= 'orange', alpha= 0.35, label= 'Aesthenosphere');
    ax.legend();
    ax.set_ylabel("Depth (km)", fontsize= 14);
    ax.set_xlabel("Temperature (K)", fontsize= 14);
    return ax;


