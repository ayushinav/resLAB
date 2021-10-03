import numpy as np;
import matplotlib.pyplot as maps;

def melt_MM3(X, P= 1):
    # fertile
    C_melt= 205.36*X*X -95.41*X+ 19.717;
    if P==2:
        C_melt= C_melt/(10**0.2);
    return C_melt;

def melt_DMM1(X, P= 1):
    # depleted
    C_melt= -5807.1*X**3+ 1495.2*X**2- 105.29*X+ 10.045;
    if P==2:
        C_melt= C_melt/(10**0.2);
    return C_melt;

def melt_T1(X, P= 1):
    # depleted
    C_melt= -2496.1*X**3+ 1320.1*X**2- 231.76*X+ 21.449;
    if P==2:
        C_melt= C_melt/(10**0.2);
    return C_melt;
    
def melt_basaltic(T, P= 1):
    R= 8.314 # J/mol/K
    C_melt= 2.321* 10**5* np.exp(-1.40* 10**5/R/T);
    if P==2:
        C_melt= C_melt/(10**0.2);
    return C_melt;

def bulk_cond_melt(X, melt_cond, solid_cond):
    HS_u= melt_cond+ (1- X)/(1/(solid_cond+ melt_cond)+
                            X/(3*melt_cond));
    HS_l= solid_cond+ (1- X)/(1/(solid_cond+ melt_cond)+
                            (1- X)/(3*melt_cond));
    
    HS= {"HS_upper": HS_u, "HS_lower": HS_l};
    return HS;
