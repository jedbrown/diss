#!/usr/bin/env python2

from __future__ import division

scalar = 8
index  = 4

class matrix_kernel:
    def bandwidth(self, rownz, flops):
        return flops / self.intensity(rownz)
    def bandwidth_fraction(self, rownz, flops, bandwidth):
        return self.bandwidth(rownz, flops) / bandwidth

class aij_no_inode(matrix_kernel):
    def intensity(self, n):          # flops / byte
        return 2*n / ((n+1)*scalar + (n+1)*index)

class aij_inode(matrix_kernel):
    def __init__(self, b): self.b = b
    def intensity(self, n):
        b = self.b
        return 2*n*b/((n+1)*b*scalar + (n+2)*index)

class baij(matrix_kernel):
    def __init__(self, b): self.b = b
    def intensity(self, n):
        b = self.b
        return 2*n*b / ((n+1)*b*scalar + (n/b+1)*index)

class sbaij(matrix_kernel):
    def __init__(self, b): self.b = b
    def intensity(self, n):
        b = self.b
        return (4*n*b - 2*b**2) / ((n+1)*b*scalar + (n/b+1)*index)

class Hardware:
    def __init__(self, bandwidth): self.bandwidth = bandwidth
    def table(self, a, b, c, d):
        bs = 2
        n = 27*bs
        bandwidth = self.bandwidth
        return ['%4d(%.0f\\%%)' % (flops,fraction*100) for (flops,fraction) in (
                (a, aij_no_inode().bandwidth_fraction(n, a, bandwidth)),
                (b, aij_inode(bs).bandwidth_fraction(n, b, bandwidth)),
                (c, baij(bs).bandwidth_fraction(n, c, bandwidth)),
                (d, sbaij(bs).bandwidth_fraction(n, d, bandwidth)),
                )]

def arithmetic_intensity_tensor(b,n):
    'n = p+1, b = block size'
    fma_lib     = 9*(3+2*b)*n**4 + (31+18*b)*n**3
    fma_physics = (3*b*3 + 4*b + b**2)*n**3     # dv : (du + dw (dw : du)) + v . w . div u + v . dw u
    fma_all     = fma_lib + fma_physics
    stored_mem  = (3+1)*b * n**3         # Store a gradient and the function value
    memory_in   = (b+3)*n**3
    unique_in   = (b+3)*(n-1)**3
    memory_out  = b*n**3
    unique_out  = b*(n-1)**3
    elem_map    = n**3
    all_mem     = scalar*(memory_in + memory_out) + index*elem_map + stored_mem*scalar
    unique_mem  = scalar*(unique_in + unique_out) + index*elem_map + stored_mem*scalar
    #my_out = memory_out; my_mem = all_mem
    my_out = unique_out; my_mem = unique_mem
    return (2*fma_all / my_out,         # flops / result
            my_mem / my_out,            # bytes / result
            2*fma_all / my_mem)             # flops / byte

def arithmetic_intensity_assembled(b,n):
    thiselem = [8 * (2*n-1)**3,                  # corners
                12 * (n-2) * n * (2*n-1)**2,     # edges
                6 * (n-2)**2 * n**2 * (2*n-1),   # faces
                (n-2)**3 * n**3]
    sharing = [8, 4, 2, 1]
    meanbrow = sum(t/s for t,s in zip(thiselem,sharing)) / (n-1)**3
    meanrow = meanbrow * b
    intensity = baij(b).intensity(meanrow)
    return (2*meanrow, 2*meanrow / intensity, intensity)


def plot_intensity(opts):
    import numpy as np
    import matplotlib.pyplot as plt
    def set_sizes_paper():
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
    set_sizes_paper()
    if opts.format != 'native': plt.rc(('backend', opts.format))
    print('tensor   :', arithmetic_intensity_tensor(1,3))
    print('assembled:', arithmetic_intensity_assembled(1,3))
    nrange = np.arange(2,9)
    brange = np.arange(1,6,2)
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1e = ax1.twinx()
    ax2e = fig.add_subplot(122)
    ax2 = ax2e.twinx()
    styles = ['-', '--', '-.']
    for i,b in enumerate(brange):
        tensor = np.array([arithmetic_intensity_tensor(b,n) for n in nrange])
        assemb = np.array([arithmetic_intensity_assembled(b,n) for n in nrange])
        ax1.semilogy(nrange-1, tensor[:,1], 's', linestyle=styles[i], label='tensor $b=%d$'%b)
        ax1.semilogy(nrange-1, assemb[:,1], 'o', linestyle=styles[i], label='assembled $b=%d$'%b)
        ax2.semilogy(nrange-1, tensor[:,0], 's', linestyle=styles[i], label='tensor $b=%d$'%b)
        ax2.semilogy(nrange-1, assemb[:,0], 'o', linestyle=styles[i], label='assembled $b=%d$'%b)
    ax1.set_xlabel('polynomial order')
    ax2.set_xlabel('polynomial order')
    ax2e.set_xlabel('polynomial order')
    ax1.set_ylabel('bytes/result')
    ax2.set_ylabel('flops/result')
    idx = [0,2,4,1,3,5]
    legend = plt.figlegend((ax1.lines[i] for i in idx),
                           (ax1.lines[i].get_label() for i in idx),
                           #(0.555,0.6), frameon=False)
                           (0.68,0.12), frameon=True)
    plt.axes(ax2); plt.xticks(nrange-1); plt.ylim(35)
    plt.axes(ax1); plt.xticks(nrange-1); plt.xlim(min(nrange-1),max(nrange-1)); plt.ylim(60,50e3)
    plt.axes(ax2e); plt.xticks(nrange-1); plt.xlim(min(nrange-1),max(nrange-1)); plt.yticks([])
    plt.axes(ax1e); plt.xticks(nrange-1); plt.xlim(min(nrange-1),max(nrange-1)); plt.yticks([])
    #plt.show()
    plt.savefig(opts.output, dpi=600)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='SpMV performance model')
    parser.add_argument('--spmv-table', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--format', choices='native png jpg eps pdf'.split(), help='Output format for plotting')
    parser.add_argument('--width-pt', help='Width, in LaTeX points, of the figure', type=float)
    parser.add_argument('-o', '--output', help='Output filename')
    opts = parser.parse_args()
    if opts.spmv_table:
        core2 = Hardware(5420).table(916, 855, 1013, 1429)
        opteron = Hardware(5787).table(506, 592, 735, 744)
        opteron4 = Hardware(14800).table(1673, 1978, 2450, 3089)
        opteron8 = Hardware(13100).table(2374, 2408, 2897, 4331)
        print('core2  :', ' & '.join(core2))
        print('opteron:', ' & '.join(opteron))
        print('opteron4:', ' & '.join(opteron4))
        print('opteron8:', ' & '.join(opteron8))
    if opts.plot:
        plot_intensity(opts)

if __name__ == "__main__":
    try: main()
    except:
        import pdb, traceback, sys
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
