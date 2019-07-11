#!/usr/bin/env python
"""
Test of the Python wrapper to the DISORT library

Module '_disort' is auto-generated with f2py (version:2).
"""

import numpy as np
import disort

##########################################################################################################

if __name__ == '__main__':

    uTau   = np.array([0.,1.])
    phi    = np.array([0.,60.,120.])
    fbeam  = 1.
    umu0   = 1./np.sqrt(2.)
    phi0   = 0.0
    albedo = 0.1
    umu    = np.array([-1.,-0.5,0.5,1.])

    # dTau   = np.ones(50)*1./50
    # w0     = np.ones(50)*1.
    # iphas  = np.ones(50,dtype='int')*2
    # gg     = np.ones(50)*0.85

    # [rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed] =\
    #                                 disort.run(dTau, w0=w0, iphas=iphas, gg=gg,
    #                                            umu0=umu0, phi0=phi0, albedo=albedo, fbeam=fbeam,
    #                                            utau=uTau, umu=umu, phi=phi)

    # CASE IMPLEMENTED IN run_disort.f
    dTau  = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    iphas = np.array([2, 2, 2, 3, 2, 2])
    w0    = np.array([0.5, 0.5, 0.5, 0.899999976, 0.5, 0.5])
    N = len(dTau)
    gg     = np.ones(N)*0.85
    prnt   = np.array([True, True, True, False, True])

    [rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed] =\
                                    disort.run(dTau, w0=w0, iphas=iphas, gg=gg,
                                               umu0=umu0, phi0=phi0, albedo=albedo, fbeam=fbeam,
                                               utau=uTau, umu=umu, phi=phi, prnt=prnt)

    print '\n\n'
    print '########################################### REFERENCE ###############################################'
    print '\n\n'
    print ' ****************************************************************************************************'
    print ' DISORT: Test Case No. 10a:  like 9c, USRANG = True                                                  '
    print ' ****************************************************************************************************'
    print ''
    print ' No. streams =  32     No. computational layers =   6'
    print '   2 User optical depths :    0.0000    1.0000'
    print '   4 User polar angle cosines : -1.00000 -0.50000  0.50000  1.00000'
    print '   3 User azimuthal angles :     0.00    60.00   120.00'
    print ' No thermal emission'
    print ' Boundary condition flag: IBCND = 0'
    print '    Incident beam with intensity =  1.000E+00 and polar angle cosine =  0.70711  and azimuth angle =   0.00'
    print '    plus isotropic incident intensity =  0.000E+00'
    print '    Bottom albedo (Lambertian) =  0.1000'
    print ' Uses delta-M method'
    print ' Uses TMS/IMS method'
    print ' Calculate fluxes and intensities'
    print ' Relative convergence criterion for azimuth series =   0.00E+00'
    print ''
    print '                                     <------------- Delta-M --------------->'
    print '                  Total    Single                           Total    Single'
    print '      Optical   Optical   Scatter   Separated   Optical   Optical   Scatter    Asymm'
    print '        Depth     Depth    Albedo    Fraction     Depth     Depth    Albedo   Factor'
    print '  1    0.0000    0.0000   0.50000     0.00000    0.0000    0.0000   0.50000   0.0000'
    print '  2    0.0000    0.0000   0.50000     0.00000    0.0000    0.0000   0.50000   0.0000'
    print '  3    0.0000    0.0000   0.50000     0.00000    0.0000    0.0000   0.50000   0.0000'
    print '  4    1.0000    1.0000   0.90000     0.00551    0.9950    0.9950   0.89950   0.8500'
    print '  5    0.0000    1.0000   0.50000     0.00000    0.0000    0.9950   0.50000   0.0000'
    print '  6    0.0000    1.0000   0.50000     0.00000    0.0000    0.9950   0.50000   0.0000'
    print ''
    print ' Number of Phase Function Moments =    33'
    print ' Layer   Phase Function Moments'
    print '    1   1.000000   0.000000   0.100000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000'
    print '    2   1.000000   0.000000   0.100000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000'
    print '    3   1.000000   0.000000   0.100000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000'
    print '    4   1.000000   0.850000   0.722500   0.614125   0.522006   0.443705   0.377150   0.320577   0.272491   0.231617'
    print '        0.196874   0.167343   0.142242   0.120906   0.102770   0.087354   0.074251   0.063113   0.053646   0.045599'
    print '        0.038760   0.032946   0.028004   0.023803   0.020233   0.017198   0.014618   0.012425   0.010562   0.008977'
    print '        0.007631   0.006486   0.005513'
    print '    5   1.000000   0.000000   0.100000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000'
    print '    6   1.000000   0.000000   0.100000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
    print '        0.000000   0.000000   0.000000'
    print '\n'
    print '                    <----------------------- FLUXES ----------------------->'
    print '  Optical  Compu    Downward    Downward    Downward      Upward                    Mean      Planck   d(Net Flux)'
    print '    Depth  Layer      Direct     Diffuse       Total     Diffuse         Net   Intensity      Source   / d(Op Dep)'
    print ''
    print '   0.0000      1   7.071E-01   1.192E-07   7.071E-01   8.424E-02   6.229E-01   9.523E-02   0.000E+00     5.983E-01'
    print '   1.0000      4   1.719E-01   3.888E-01   5.607E-01   5.607E-02   5.046E-01   8.305E-02   0.000E+00     1.044E-01'
    print '\n'
    print ' *********  I N T E N S I T I E S  *********'
    print ''
    print '             Polar   Azimuth angles (degrees)'
    print '  Optical   Angle'
    print '   Depth   Cosine'
    print ''
    print '                         0.00      60.00     120.00'
    print '   0.0000 -1.0000  0.000E+00  0.000E+00  0.000E+00'
    print '          -0.5000  0.000E+00  0.000E+00  0.000E+00'
    print '           0.5000  4.465E-02  3.413E-02  2.211E-02'
    print '           1.0000  2.149E-02  2.149E-02  2.149E-02'
    print ''
    print '                         0.00      60.00     120.00'
    print '   1.0000 -1.0000  4.352E-02  4.352E-02  4.352E-02'
    print '          -0.5000  6.599E-01  7.163E-02  2.009E-02'
    print '           0.5000  1.785E-02  1.785E-02  1.785E-02'
    print '           1.0000  1.785E-02  1.785E-02  1.785E-02'
