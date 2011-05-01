#!/usr/bin/env python2

from collections import namedtuple
from pylab import *

Sequence = namedtuple('Sequence', 'order samples')
Sample = namedtuple('Sample', 'nelem l2 lmax h1 hmax')

q1 = Sequence(order=1, samples=[Sample(1, 1.79e+00, 6.50e-01, 3.70e+00, 1.08e+00),
                                Sample(2, 5.49e-01, 3.40e-01, 1.61e+00, 6.92e-01),
                                Sample(4, 1.53e-01, 1.26e-01, 8.01e-01, 4.51e-01),
                                Sample(8, 3.94e-02, 3.73e-02, 3.98e-01, 2.81e-01),
                                Sample(16, 9.95e-03, 1.01e-02, 1.98e-01, 1.57e-01),
                                Sample(32, 2.49e-03, 2.61e-03, 9.92e-02, 8.32e-02),
                                ])
q2 = Sequence(order=2, samples=[Sample(1, 2.44e-01, 1.82e-01, 9.48e-01, 4.60e-01),
                                Sample(2, 3.71e-02, 4.47e-02, 2.86e-01, 1.54e-01),
                                Sample(4, 4.48e-03, 6.23e-03, 6.94e-02, 4.34e-02),
                                Sample(8, 5.60e-04, 9.31e-04, 1.74e-02, 1.29e-02),
                                Sample(16, 7.01e-05, 1.23e-04, 4.34e-03, 3.52e-03),
                                ])
q3 = Sequence(order=3, samples=[Sample(1, 4.14e-02, 2.71e-02, 2.90e-01, 1.63e-01),
                                Sample(2, 2.06e-03, 2.06e-03, 2.39e-02, 1.14e-02),
                                Sample(4, 1.81e-04, 2.06e-04, 4.23e-03, 2.88e-03),
                                Sample(8, 1.22e-05, 1.87e-05, 5.79e-04, 5.84e-04),
                                ])
q5 = Sequence(order=5, samples=[Sample(1, 3.76e-03, 2.90e-03, 4.69e-02, 3.16e-02),
                                Sample(2, 7.58e-05, 5.92e-05, 1.62e-03, 1.05e-03),
                                Sample(4, 7.33e-07, 6.61e-07, 2.59e-05, 1.76e-05),
                                #Sample(8, 2.48e-08, 1.06e-07, 1.84e-06, 1.13e-05), # lots of quadrature error
                                ])
q7 = Sequence(order=7, samples=[Sample(1, 4.46e-04, 3.59e-04, 8.15e-03, 5.83e-03),
                                Sample(2, 2.95e-06, 2.95e-06, 8.21e-05, 6.05e-05),
                                Sample(4, 7.65e-09, 1.07e-08, 4.09e-07, 3.95e-07),
                                ])
q9 = Sequence(order=9, samples=[Sample(1, 5.81e-05, 5.04e-05, 1.42e-03, 1.05e-03),
                                Sample(2, 6.27e-08, 7.59e-08, 1.63e-06, 1.60e-06)
                                ])

class Container:
    pass

def table(seq):
    for i,s in enumerate(seq.samples):
        o = Container()
        if i == 0:
            for attr in 'l2 lmax h1 hmax'.split():
                setattr(o, attr, '---')
        else:
            def rate(a, b):
                return '%3.2f' % (-log2(b/a),)
            sm1 = seq.samples[i-1]
            for attr in 'l2 lmax h1 hmax'.split():
                setattr(o,attr, rate(getattr(sm1,attr), getattr(s,attr)))
        yield '$Q_%d$ & $%d^3$ & %d & %5.2e & %s & %5.2e & %s & %5.2e & %s & %5.2e & %s' % (seq.order, s.nelem, (s.nelem*seq.order+1)**3, s.l2, o.l2, s.lmax, o.lmax, s.h1, o.h1, s.hmax, o.hmax)

print (r'\\' + '\n' + r'\hline' + '\n').join([(r' \\'+ '\n').join(table(q)) for q in [q1, q2, q3, q5, q7, q9]])
