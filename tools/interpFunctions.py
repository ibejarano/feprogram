import numpy as np

def funH4(r,s):
    r_plus = 1+r
    r_minus = 1-r
    s_plus = 1+s
    s_minus = 1-s

    h4 = 0.25*r_plus*s_minus
    h3 = 0.25*r_minus*s_minus
    h2 = 0.25*r_minus*s_plus
    h1 = 0.25*r_plus*s_plus

    H=np.matrix([h1,h2,h3,h4])
    return H

def funH9(r,s):
    r_power = 1-r**2
    s_power = 1-s**2
    r_plus = 1+r
    r_minus = 1-r
    s_plus = 1+s
    s_minus = 1-s

    h9 = r_power*s_power
    h8 = 0.5*s_power*r_plus - 0.5*h9
    h7 = 0.5*r_power*s_minus - 0.5*h9
    h6 = 0.5*s_power*r_minus - 0.5*h9
    h5 = 0.5*r_power*s_plus - 0.5*h9
    h4 = 0.25*r_plus*s_minus - 0.5*h7 - 0.5*h8 -0.25*h9
    h3 = 0.25*r_minus*s_minus - 0.5*h6 - 0.5*h7 - 0.25*h9
    h2 = 0.25*r_minus*s_plus - 0.5*h5 - 0.5*h6 - 0.25*h9
    h1 = 0.25*r_plus*s_plus - 0.5*h5 -0.5*h8 -0.25*h9

    H=np.matrix([h1,h2,h3,h4,h5,h6,h7,h8,h9])

    return H

def funHrs4(r,s):
    hrquad4 = [1+s , -1-s , -1+s ,1-s]
    hsquad4 = [1+r,1-r , -1 +r , -1 -r]

    hr = [0]*4
    hs = [0]*4
    for i in range(4):
        hr[i] = 0.25*hrquad4[i]
        hs[i] = 0.25*hsquad4[i]

    dhrs = np.matrix(np.zeros((2,4)))

    dhrs[0] = hr
    dhrs[1] = hs
    return dhrs    

def funHrs9(r,s):
    hrquad4 = [1+s , -1-s , -1 +s ,1-s]
    hrquad8 = [-2*r*(1+s),-1*(1-s**2),-2*r*(1-s),(1-s**2)]
    hrquad9 = -2*r*(1-s**2)
    hsquad4 = [1+r,1-r ,-1+r ,-1-r]
    hsquad8 = [(1-r**2),-2*s*(1-r),-1*(1-r**2),-2*s*(1+r)]
    hsquad9 = -2*s*(1-r**2)

    hr = [0]*9
    hs = [0]*9

    for i in range(4):
        hr[i+4] = 0.5*hrquad8[i] - 0.5*hrquad9
        hs[i+4] = 0.5*hsquad8[i] - 0.5*hsquad9

    hr[8] = hrquad9
    hs[8] = hsquad9
    
    for i in range(1,4):
        hr[i] = 0.25*hrquad4[i] - 0.5*hr[i+4] - 0.5*hr[i+3] - 0.25*hrquad9
        hs[i] = 0.25*hsquad4[i] - 0.5*hs[i+4] - 0.5*hs[i+3] - 0.25*hsquad9
        
    hr[0] = 0.25*hrquad4[0] - 0.5*hr[4] - 0.5*hr[7] - 0.25*hrquad9
    hs[0] = 0.25*hsquad4[0] - 0.5*hs[4] - 0.5*hs[7] - 0.25*hsquad9
    
    dhrs = np.matrix(np.zeros((2,9)))

    dhrs[0] = hr
    dhrs[1] = hs
    return dhrs    
