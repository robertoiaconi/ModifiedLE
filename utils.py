import numpy as np
import matplotlib.pyplot as plt

# Initialise physical constants
g_Msun  = 1.989e33
cm_Rsun = 6.957e10
cgs_G   = 6.6743e-8
cgs_kb  = 1.38e-16
cgs_amu = 1.66e-24 


def smooth(x,window_len=11,window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]

    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    return y[(window_len//2-1):-(window_len//2)]

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

class Kernel:
    def __init__(self, name, radius, func):
        self.name = name
        self.radius = radius
        self.func   = func

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def set_hsoft(self, h):
        self.h = h

    def kernel(self, r, h=None):
        if h is not None:
            return self.func(r,h)
        else:
            return self.func(r,self.h)

    # This could do with some prettifying
    def plot_kernel(self):
        fig, axes = plt.subplots(1,2)

        rs      = np.linspace(0, self.radius, 1000)
        us      = np.array([])
        chis    = np.array([])
        dchidus = np.array([])

        for r in rs:
            u, chi, dchidu = self.kernel(r)
            us = np.append(us,u)
            chis = np.append(chis,chi)
            dchidus = np.append(dchidus,dchidu)

        axes[0].plot(rs, us*chis)
        axes[1].plot(rs, chis + us / 3. * dchidus)
        plt.show()


def ohlmann_kernel(r, h):
    u = r / h
    if (u >= 0 and u < 0.5):
        chi    = 32. / 3. - 192. * u**2 / 5. + 32. * u**3
        dchidu = -384.*u / 5. + 96. * u**2
    elif (u >= 0.5 and u < 1.): 
        chi    = -1. / (15. * u**3) + 64. / 3. - 48. * u + 192. * u**2 / 5. - 32. * u**3 / 3
        dchidu = 1. / (5 * u**4) - 48. + 384. * u / 5 - 32. * u**2
    else:
        chi    = 1. / u**3
        dchidu = -3. / u**4

    return u, chi, dchidu

def phantom_kernel(r, h):
    u = r / h

    if (u >= 0 and u < 1.):
        chi    = (15.*u**3 - 36.*u**2 + 40.)/30.
        dchidu = (15.*u**2 - 24.*u)/10.
    elif (u >= 1. and u < 2.):
        chi    = (-5.*u**5 + 36.*u**4 - 90.*u**3 + 80.*u**2 - 2.*u**(-1))/(30.*u**2)
        dchidu = (-5.*u**6 + 24.*u**5 - 30.*u**4 + 2)/(10.*u**4)
    else:
        chi    = 1./u**3
        dchidu = -3./u**4
    return u, chi, dchidu

ohlmann = Kernel('ohlmann', 1, ohlmann_kernel)
phantom = Kernel('phantom', 2, phantom_kernel)

