# Importing necessary libraries

import numpy as np;
import matplotlib.pyplot as maps;
import math;

def ol_cond(T, Cw=0.0):
    """
    T= Temperature in Kelvin
    Cw= Concentration of water
    """
    k= 8.61734279*10**-5;
    # ionic 
    A_ionic= 10**4.73;
    H_ionic= 2.31;
    
    # polaron
    A_polaron= 10**2.98;
    H_polaron= 1.71;
    
    # proton
    A_proton= 10**1.90;
    H_proton= 0.92;
    alpha= 0.16;
    
    idx_ionic= (T>= 1700);
    idx_polaron= np.logical_and(T>=1300,T<1700);
    idx_proton= (T< 1300);
    
    C_ionic= A_ionic* np.exp(-H_ionic/k/T);
    C_polaron= A_polaron* np.exp(-H_polaron/k/T);
    C_proton= A_proton* Cw* np.exp(-(H_proton- alpha* Cw**(1/3))/k/T);
    
    C_tot= C_ionic+ C_polaron+ C_proton;
    
    id_ionic= np.array(idx_ionic, dtype= int);
    id_polaron= np.array(idx_polaron, dtype= int);
    id_proton= np.array(idx_proton, dtype= int);
    
    C_ol= C_ionic*id_ionic+ C_polaron*id_polaron+ C_proton*id_proton;
    
    
    idx= {'ionic': idx_ionic, 'polaron': idx_polaron, 'proton': idx_proton};
    C= {'ionic': C_ionic, 'polaron': C_polaron, 'proton': C_proton, 'total': C_tot, 'ol': C_ol};
    
    return [idx, C];


def opx_cond(T, Cw=0.0):
    """
    T= Temperature in Kelvin
    Cw= Concentration of water
    """
    k= 8.61734279*10**-5;
    
    # polaron
    A_polaron= 10**3.99;
    H_polaron= 1.88;
    
    # proton
    A_proton= 10**2.58;
    H_proton= 0.84;
    alpha= 0.08;
    
    idx_polaron= (T>=1300);
    idx_proton= (T< 1300);
    
    C_polaron= A_polaron* np.exp(-H_polaron/k/T);
    C_proton= A_proton* Cw* np.exp(-(H_proton- alpha* Cw**(1/3))/k/T);
    
    C_tot= C_polaron+ C_proton;
    
    id_polaron= np.array(idx_polaron, dtype= int);
    id_proton= np.array(idx_proton, dtype= int);
    
    C_opx= C_polaron*id_polaron+ C_proton*id_proton;
    
    idx= {'polaron': idx_polaron, 'proton': idx_proton};
    C= {'polaron': C_polaron, 'proton': C_proton, 'total': C_tot, 'opx': C_opx};
    
    return [idx, C];


def cpx_cond(T, Cw=0.0):
    """
    T= Temperature in Kelvin
    Cw= Concentration of water
    """
    k= 8.61734279*10**-5;
    
    # polaron
    A_polaron= 10**2.16;
    H_polaron= 1.06 # 102 kJ/mol
    
    # proton
    A_proton= 10**3.56;
    H_proton= 0.73; # 71 kJ/mol
    r= 1.13;
    
    idx_polaron= (T>=0);
    idx_proton= (T>=0);
    
    C_polaron= A_polaron* np.exp(-H_polaron/k/T);
    C_proton= A_proton* (Cw**r)* np.exp(-H_proton/k/T);
    
    C_tot= C_polaron+ C_proton;
    
    id_polaron= np.array(idx_polaron, dtype= int);
    id_proton= np.array(idx_proton, dtype= int);
    if Cw==0:
        C_cpx= C_polaron;
    else:
        C_cpx= C_proton;
        
    idx= {'polaron': idx_polaron, 'proton': idx_proton};
    C= {'polaron': C_polaron, 'proton': C_proton, 'total': C_tot, 'cpx': C_cpx};
    
    return [idx, C];


def gt_cond(T, Cw=0.0):
    """
    T= Temperature in Kelvin
    Cw= Concentration of water
    """
    k= 8.61734279*10**-5;
    P= 8;
    # polaron
    A0= 1036;
    B= 0.044;
    A_polaron= A0*(1-B*P);
    E_polaron= 128;
    V_polaron= 2.50;
    H_polaron= (E_polaron+ P*V_polaron)/96.48 # kJ/mol => eV
    
    # proton
    A_proton= 1950;
    E_proton= 70;
    V_proton= -0.57;
    H_proton= (E_proton+ P*V_proton)/96.48; # kJ/mol => eV
    r= 0.63;
    
    idx_polaron= (T>=0);
    idx_proton= (T>=0);
    
    C_polaron= A_polaron* np.exp(-H_polaron/k/T);
    C_proton= A_proton* (Cw**r)* np.exp(-H_proton/k/T);
    
    C_tot= C_polaron+ C_proton;
    
    id_polaron= np.array(idx_polaron, dtype= int);
    id_proton= np.array(idx_proton, dtype= int);
    if Cw==0:
        C_gt= C_polaron;
    else:
        C_gt= C_proton;
        
    idx= {'polaron': idx_polaron, 'proton': idx_proton};
    C= {'polaron': C_polaron, 'proton': C_proton, 'total': C_tot, 'gt': C_gt};
    
    return [idx, C];
