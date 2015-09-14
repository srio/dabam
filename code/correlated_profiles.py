"""

aspheric_fit
   tests for rough surface generation & analysis

   inputs: P

    output: some plots

    Functions translated from:

       http://www.mysimlabs.com/surface_generation.html

 
       MODIFICATION HISTORY:
           20150910 srio@esrf.eu, written
"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


import numpy
from matplotlib import pylab as plt



def rsgeng1D(N,rL,h,cl):
    # function [f,x] = rsgeng1D(N,rL,h,cl)
    # %
    # % [f,x] = rsgeng1D(N,rL,h,cl)
    # %
    # % generates a 1-dimensional random rough surface f(x) with N surface points.
    # % The surface has a Gaussian height distribution function and a Gaussian
    # % autocovariance function, where rL is the length of the surface, h is the
    # % RMS height and cl is the correlation length.
    # %
    # % Input:    N   - number of surface points
    # %           rL  - length of surface
    # %           h   - rms height
    # %           cl  - correlation length
    # %
    # % Output:   f   - surface heights
    # %           x   - surface points
    # %
    # % Last updated: 2010-07-26 (David Bergström).
    # %
    #

    # format long;
    #
    # x = linspace(-rL/2,rL/2,N);
    x = numpy.linspace(-rL/2,rL/2,N)

    #
    # Z = h.*randn(1,N); % uncorrelated Gaussian random rough surface distribution
    #                      % with mean 0 and standard deviation h
    #
    Z = h * numpy.random.randn(1.0,N)
    Z.shape = -1

    # % Gaussian filter
    # F = exp(-x.^2/(cl^2/2));
    F = numpy.exp(-x**2/(cl**2/2))
    #
    # % correlation of surface using convolution (faltung), inverse
    # % Fourier transform and normalizing prefactors
    # f = sqrt(2/sqrt(pi))*sqrt(rL/N/cl)*ifft(fft(Z).*fft(F));
    f = numpy.sqrt(2/numpy.sqrt(numpy.pi))*numpy.sqrt(rL/N/cl)*numpy.fft.ifft(numpy.fft.fft(Z)*numpy.fft.fft(F))

    return x,f


def acf1D(f,x):
    # function [acf,cl,lags] = acf1D(f,x,opt)
    # %
    # % [acf,cl,lags] = acf1D(f,x)
    # %
    # % calculates the autocovariance function and correlation length of a
    # % 1-d surface profile f(x).
    # %
    # % Input:    f    - surface heights
    # %           x    - surface points
    # %           opt - optional parameter (type 'plot' for plotting the
    # %                 normalized autocovariance function)
    # %
    # % Output:   acf  - autocovariance function
    # %           cl   - correlation length
    # %           lags - lag length vector (useful for plotting the acf)
    # %
    # % Last updated: 2010-07-26 (David Bergstrom)
    # %

    #
    # format long
    #
    # N = length(x); % number of sample points
    N = len(x)

    # lags = linspace(0,x(N)-x(1),N); % lag lengths
    lags = numpy.linspace(0,x[-1]-x[0],N)

    #
    # % autocovariance function calculation
    # c=xcov(f,'coeff'); % the autocovariance function
    f -= f.mean()
    c = numpy.convolve(f,f[::-1])
    c = c / c.max()

    # acf=c(N:2*N-1); % right-sided version
    acf=c[(N-1):2*N-2]

    #
    # % correlation length calculation
    # k = 1;
    k = 0

    while acf[k] > 1/numpy.exp(1):
        k = k + 1

    # while (acf(k) > 1/exp(1))
    #     k = k + 1;
    # end;
    # cl = 1/2*(x(k-1)+x(k)-2*x(1)); % the correlation length
    #
    cl = 1/2*(x[k-1]+x[k]-2*x[0])

    # % optional plotting
    # if nargin<3 || isempty(opt)
    #     return;
    # end;
    # if nargin==3
    #     if ischar(opt)
    #         if strcmp(opt,'plot');
    #             plot(lags,acf);
    #             xlabel('lag length')
    #             ylabel('acf (normalized)')
    #             title('Plot of the normalized acf')
    #         else fprintf('%s is not a valid option. Type \''help acf1D\'' for further details.\n',opt);
    #         end;
    #     else fprintf('Option must be a string. Type \''help acf1D\'' for further details.\n');
    #     end;
    # end;
    return acf,cl,lags

def psd(x, y, onlyrange = None):
    """
     psd: Calculates the PSD (power spectral density) from a profile

      INPUTS:
           x - 1D array of (equally-spaced) lengths.
           y - 1D array of heights.
      OUTPUTS:
           f - 1D array of spatial frequencies, in units of 1/[x].
           s - 1D array of PSD values, in units of [y]^3.
      KEYWORD PARAMETERS:
           onlyrange - 2-element array specifying the min and max spatial
               frequencies to be considered. Default is from
               1/(length) to 1/(2*interval) (i.e., the Nyquist
               frequency), where length is the length of the scan,
               and interval is the spacing between points.

      PROCEDURE
            S=Length*ABS(ifft(Y*Window)^2
            Where Length is as described above, and Window is the value of
            the optional window function

    """
    n_pts = x.size
    if (n_pts <= 1):
        print ("psd: Error, must have at least 2 points.")
        return 0

    xx = x
    yy = y

    window=yy*0+1.
    length=xx.max()-xx.min()  # total scan length.
    # psd
    s=length*numpy.absolute(numpy.fft.ifft(yy*window)**2)

    s=s[1:(n_pts/2+1*numpy.mod(n_pts,2))]  # take an odd number of points.
    n_ps=s.size                       # number of psd points.
    interval=length/(n_pts-1)         # sampling interval.
    f_min=1./length                   # minimum spatial frequency.
    f_max=1./(2.*interval)            # maximum (Nyquist) spatial frequency.
    # spatial frequencies.
    f=numpy.arange(float(n_ps))/(n_ps-1)*(f_max-f_min)+f_min

    if onlyrange != None:
        roi =  (f <= onlyrange[1]) * (f >= onlyrange[0])
        if roi.sum() > 0:
            roi = roi.nonzero()
            f = f[roi]
            s = s[roi]

    return s,f


def main():
    #

    #
    # create profile
    #

    N = 1000
    mirror_length = 1.0
    height_rms = 3e-9

    profile_type = 0 # 0=Fractal, 1=Gaussian

    if profile_type == 0:
        xmirr = numpy.linspace(-0.5*mirror_length,0.5*mirror_length,N)
        xstep = mirror_length/(N-1)

        freq = numpy.linspace(1/(1*mirror_length),1/(4*xstep),500)
        print("frequency range from %f, to %f"%(freq.min(),freq.max()))
        ampl = freq**(-0.9)
        phases = numpy.random.rand(freq.size)*2*numpy.pi
        ymirr = numpy.zeros(N)
        for i in range(len(freq)):
            ymirr += (ampl[i] *  numpy.sin(2*numpy.pi*freq[i]*xmirr + phases[i]))

        x = xmirr
        f = ymirr / ymirr.std() * height_rms
    elif profile_type == 1:
        height_rms = 3e-9
        correlation_length = 0.03
        x, f = rsgeng1D(N,mirror_length,height_rms,correlation_length)

    acf,cl,lags = acf1D(f,x)
    #print("acf=",acf,"cl=",cl,"lags=",lags)
    if profile_type == 0:
        print("f RMS in nm: obtained:%f"%(1e9*f.std()))
        print("f cl in m: obtained:%f"%(cl))
    else:
        print("f RMS in nm: defined=%f, obtained:%f"%(1e9*height_rms,1e9*f.std()))
        print("f cl in m: defined=%f, obtained:%f"%(correlation_length,cl))

    psd_y, psd_x = psd(x,f)


    outFile = ""
    if outFile != "":
        dd=numpy.concatenate( (x, f) ,axis=0).reshape(2,-1).transpose()
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk:\n")

    #
    #plots
    #
    do_plots = 1
    if do_plots:
        #
        #plots
        #
        from matplotlib import pylab as plt

        plt.figure(1)
        plt.plot(x,f)
        plt.title("Simulated profile")
        plt.xlabel("Y [m]")
        plt.ylabel("Z [m]")
        plt.show()

        plt.plot(lags[0:-1],acf)
        plt.title("Autocorrelation function")
        plt.xlabel("Length [m]")
        plt.ylabel("autocorrelation")
        plt.show()


        plt.loglog(psd_x,psd_y)
        plt.title("PSD function")
        plt.xlabel("Spatiam frequency [m^-1]")
        plt.ylabel("PSD [m^3]")
        plt.show()

#
# main program
#
if __name__ == '__main__':
    main()
