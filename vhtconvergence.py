#!/usr/bin/env python2

from __future__ import division, with_statement
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt

def set_sizes_paper(opts):
    def get_dims(args, dflt_width, ratio=(np.sqrt(5)-1)/2):
        inches_per_pt = 1/72.27
        width = inches_per_pt * (args.width_pt if args.width_pt else dflt_width)
        height = width * ratio * 0.9
        return width, height
    fig_size = get_dims(opts, 480)
    plt.rcParams.update({'axes.titlesize': 11,
                     'axes.labelsize': 9,
                     'text.fontsize': 9,
                     'lines.linewidth': 1,
                     'legend.fontsize': 9,
                     'legend.markerscale' : 1,
                     'xtick.labelsize': 8,
                     'ytick.labelsize': 8,
                     'text.usetex': True,
                     'figure.figsize': fig_size})
    #plt.subplots_adjust(left=0.1,right=0.9,bottom=0.10,top=0.93)

Residual = namedtuple('Residual', 'mesh step all momentum mass energy')
# g 'SNES Function' vht-q323b12.log
logs = [
    Residual(12, 0 , 3.093796e+01 ,  1.431024e+01 , 2.142762e-01 , 2.742861e+01),
    Residual(12, 1 , 1.187863e+01 ,  8.386431e+00 , 1.086178e-05 , 8.412471e+00),
    Residual(12, 2 , 5.445497e+00 ,  4.103421e+00 , 2.928430e-06 , 3.579857e+00),
    Residual(12, 3 , 3.135764e+00 ,  3.059956e+00 , 1.744093e-06 , 6.853340e-01),
    Residual(12, 4 , 1.914854e+00 ,  1.891518e+00 , 8.688964e-07 , 2.980380e-01),
    Residual(12, 5 , 6.890634e-01 ,  6.852763e-01 , 4.952597e-07 , 7.214378e-02),
    Residual(12, 6 , 4.869877e-02 ,  4.827890e-02 , 1.430063e-07 , 6.381075e-03),
    Residual(12, 7 , 4.433521e-05 ,  3.257086e-05 , 1.057706e-08 , 3.007905e-05),
    Residual(12, 8 , 4.473666e-10 ,  4.249319e-10 , 1.759267e-11 , 1.387815e-10),
    ]

# for n in 2 4 6 8 12; do echo $n; awk '/nodes/{print $3, $5, $7;}' vht-q323b$n-x0diff.log; awk '/^Integral/{print $8;}' vht-q323b$n-x0diff.log; done
# mesh nu np ne u0 u1 p0 p1 e0 e1
discerrstr = """
2   512 216 512 4.85e-03 5.43e-02 1.34e-01 1.29e+00 6.57e-02 7.97e-01
4   4096 1728 4096 4.90e-04 7.27e-03 2.20e-02 2.52e-01 6.63e-03 1.60e-01
6   13824 5832 13824 4.09e-04 3.16e-03 1.81e-02 1.22e-01 1.56e-03 5.00e-02
8   32768 13824 32768 4.05e-04 2.56e-03 1.78e-02 8.95e-02 8.40e-04 2.14e-02
12  110592 46656 110592 4.05e-04 2.42e-03 1.77e-02 7.48e-02 6.96e-04 6.72e-03
16  262144 110592 262144 4.05e-04 2.41e-03 1.77e-02 7.21e-02 6.81e-04 3.46e-03
""".split('\n')[1:-1]
discerr = np.array(([map(float,line.split()) for line in discerrstr]))

def fit_ideal_loglog(x,y,slope):
    xmin, ymin = log(x[0]), log(y[0])
    b = ymin - xmin*slope;
    ymax = log(y[-1])
    xmax = (ymax - b)/slope     # Find the x associated with the rightmost data point
    fit_x = log(x) #array([xmin, xmax])
    fit_y = fit_x*slope + b
    return (exp(fit_x), exp(fit_y))
def refslope(h, order, ref):
    y = np.log(h)*order + np.log(ref)
    return np.exp(y)
def refline(plt, h, order, ref, style='k-'):
    y = refslope(h,order,ref)
    plt.loglog(h, y, style)
    plt.annotate('slope=%d'%order, (h[1],y[1]), xytext=(-20,10), textcoords='offset points')

def plot_error(opts):
    set_sizes_paper(opts)
    if opts.format != 'native': plt.rc(('backend', opts.format))
    arr = discerr
    h = 2/arr[:,0]
    plt.loglog(h, arr[:,4], 'bo-', label  = 'Momentum $L^2$')
    plt.loglog(h, arr[:,5], 'bo-.', label = 'Momentum $H^1$')
    plt.loglog(h, arr[:,6], 'gs-', label  = 'Pressure $L^2$')
    plt.loglog(h, arr[:,7], 'gs-.', label = 'Pressure $H^1$')
    plt.loglog(h, arr[:,8], 'rv-', label  = 'Energy $L^2$')
    plt.loglog(h, arr[:,9], 'rv-.', label = 'Energy $H^1$')
    refline(plt,h[:3],2, 2, 'k-')
    refline(plt,h[:3],3, 0.3, 'k-')
    refline(plt,h[:3],4, 0.02, 'k-')
    plt.legend(loc='upper left')
    plt.xlabel('Mesh size')
    plt.ylabel('Continuous norm of error')
    if opts.output: plt.savefig(opts.output, dpi=600)
    else:           plt.show()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='SpMV performance model')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--format', choices='native png jpg eps pdf'.split(), help='Output format for plotting')
    parser.add_argument('--width-pt', help='Width, in LaTeX points, of the figure', type=float)
    parser.add_argument('-o', '--output', help='Output filename')
    opts = parser.parse_args()
    if opts.plot:
        plot_error(opts)

if __name__ == "__main__":
    try: main()
    except SystemExit:
        pass
    except:
        import pdb, traceback, sys
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
