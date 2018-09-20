import numpy as np

def smooth(x,window_len=11,window='hanning'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]

    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    return y[(window_len/2-1):-(window_len/2)]

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

class Kernel:
    def __init__(self, radius, func):
        self.radius = radius
        self.func   = func

    def set_hsoft(self, h):
        self.h = h

    def kernel(self, r, h=None):
        if h is not None:
            return self.func(r,h)
        else:
            return self.func(r,self.h)

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

ohlmann = Kernel(1, ohlmann_kernel)
phantom = Kernel(2, phantom_kernel)

