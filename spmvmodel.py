#!/usr/bin/env python3

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

if __name__ == "__main__":
    core2 = Hardware(5420).table(916, 855, 1013, 1429)
    opteron = Hardware(5787).table(506, 592, 735, 744)
    print('core2  :', ' & '.join(core2))
    print('opteron:', ' & '.join(opteron))
