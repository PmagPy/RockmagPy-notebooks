#-------------------------------------------------------------------------------
# Name:        module2
# Purpose:
#
# Author:      peat
#
# Created:     23/01/2024
# Copyright:   (c) peat 2024
# Licence:     <your licence>
#-------------------------------------------------------------------------------

def main():
    pass

if __name__ == '__main__':
    main()

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as sc
import pandas as pd
from scipy.optimize import brent
from numpy.linalg import svd, lstsq, solve


def y_polate(n, xarr, yarr, xval):
    # linear interpolation/extrapolation, returns best-fit y for a specified x (e.g. for Mr)
    intercept, slope = linefit1(n, xarr, yarr)
    ycal = intercept + slope * xval

    return ycal

def x_polate(n, xarr, yarr, yval):
    intercept, slope = linefit1(n, xarr, yarr)
    xcal = (yval-intercept) / slope

    return xcal

def linefit1(n, xarr, yarr):
    sum_x = sum(xarr)
    sum_y = sum(yarr)
    sum_xy = sum(x * y for x, y in zip(xarr, yarr))
    sum_x_squared = sum(x ** 2 for x in xarr)

    slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x ** 2)
    intercept = (sum_y - slope * sum_x) / n

    return intercept, slope

def linefit2(n, xarr, yarr):
    sumx = sumy = st2 = b = ssr = sst = 0
    for i in range(n):
        sumx += xarr[i]
        sumy += yarr[i]

    xavg = sumx / n
    yavg = sumy / n

    for i in range(n):
        xdelta = xarr[i] - xavg
        st2 += xdelta ** 2
        ydelta = yarr[i] - yavg
        sst += ydelta ** 2
        b += xdelta * yarr[i]

    if st2 > 0:
        b /= st2
    else:
        b = 9E9  # slope

    a = (sumy - sumx * b) / n  # intercept

    for i in range(n):
        # ssr += (yarr[i] - (a + b * xarr[i])) ** 2  # sum squared residuals
        ssr += ((a + b * xarr[i]) - yavg) ** 2  # sum squared due to regression

    if sst > 0:
        rsqr = ssr / sst
    else:
        rsqr = 1

    slope = b
    intercept = a

    return intercept, slope, rsqr
"""  using numpy matmul instead
def mmult(a, b, nra, nca, nrb, ncb):
    if not nca == nrb:
        # print("Error: Inner matrix dimensions must agree.")
        pass
    else:
        c = [[0 for _ in range(ncb + 1)] for _ in range(nra + 1)]
        for i in range(1, nra + 1):
            for j in range(1, ncb + 1):
                c[i][j] = 0
        for i in range(1, nra + 1):
            for j in range(1, ncb + 1):
                for k in range(1, nca + 1):
                    c[i][j] += a[i][k] * b[k][j]
        return c
"""

"""  Not being used
def polynomialfit1(n, mindeg, maxdeg, xarr, yarr):  #polynomialfit1(nsmooth, 0, polydegree, x, y, y, r2)
    # y=a0x^0+a1x^1+...aNx^N
    npoints = n
    n1 = maxdeg - mindeg + 1
    npar = n1

    a = np.zeros((npoints, npar))
    SV_A = np.zeros((npoints, npar))
    sv_b = np.zeros(npoints)
    sv_x = np.zeros(npar)

    for i in range(npoints):
        for j in range(npar):
            a[i, j] = xarr[i] ** (mindeg + j)

    SV_A = a.copy()
    sv_b = yarr.copy()

    _, sv_w, SV_Vt = svd(SV_A, full_matrices=False)
    sv_w_inv = 1.0 / sv_w
    sv_x = np.dot(SV_Vt.T, np.dot(np.diag(sv_w_inv), np.dot(SV_Vt, sv_b)))

    coeffarr = sv_x.copy()

    y = np.zeros((npoints, 2))
    y[:, 1] = coeffarr

    y1 = np.dot(a, y[:, 1])

    y11 = y1.copy()

    # Assuming linefit1 is a separate function
    intercept, slope, rsqr = linefit1(npoints, yarr, y11)

    return coeffarr, rsqr
"""
def polynomialfit2(n, mindeg, maxdeg, xarr, yarr, xgiven):
    # Define coefficients array
    n1 = maxdeg - mindeg + 1
    #coeffs = np.empty(n1 + 1)
    #coeffs = np.empty(2)

    # Call polynomialfit1 to get coefficients
    coeffs = np.polyfit(xarr, yarr, maxdeg)
    #coeffs, r2 = polynomialfit1(n, mindeg, maxdeg, xarr, yarr)

    # Calculate yfit using the obtained coefficients
    yfit = 0
    j = -1   #j = 0
    for i in range(mindeg, maxdeg + 1):
        j += 1
        x = xgiven ** i
        yfit += coeffs[j] * x

    return yfit

def branch_sub(n_loop, grid_moments):

    j = n_loop // 2
    moment_sub = np.ndarray([j])
    for i in range (j):
        moment_sub[i] = grid_moments[i]-grid_moments[n_loop-1-i]
    return moment_sub


def loop_grid(n_loop, polydegree, nsmooth, loop_fields, loop_moments, B_offset, M_offset):

    n2 = nsmooth // 2  # make nsmooth odd; n2=# on each side, + 1 in center
    xvals = np.ndarray([nsmooth])
    yvals = np.ndarray([nsmooth])
    xvals[:]=0
    yvals[:]=0

    grid_fields = np.ndarray([n_loop])
    grid_moments = np.ndarray([n_loop])
    grid_fields[:] = 0
    grid_moments[:] = 0


    field_step = loop_fields[0] / n_loop * 4
    grid_fields[0] = round(1000 * (loop_fields[0] - abs(B_offset))) / 1000  # nearest mT
    grid_fields[n_loop // 2] = -round(1000 * (loop_fields[0] - abs(B_offset))) / 1000  # nearest mT
    field_step = grid_fields[0] / (n_loop - 2) * 4
    for i in range(1, n_loop // 2):  #was +1, but this is 0 bassed array
        grid_fields[i] = grid_fields[i - 1] - field_step
        grid_fields[i + n_loop // 2] = grid_fields[i - 1 + n_loop // 2] + field_step


    j = 1
    grid_moments[0] = loop_moments[0] - M_offset

    for i in range(1, n_loop // 2):
        while not (loop_fields[j] - B_offset >= grid_fields[i]): # or (j <= 2):
            if j<=1:
                break
            if j > 1:
                j -= 1
        while not (loop_fields[j] - B_offset <= grid_fields[i]): #or (j >= n_loop // 2):
            if j >=n_loop // 2:
                break
            if j < n_loop // 2:
                j += 1

        k = 1
        xvals[1] = loop_fields[j] - B_offset
        xvals[2] = loop_fields[j - k] - B_offset
        yvals[1] = loop_moments[j] - M_offset
        yvals[2] = loop_moments[j - k] - M_offset

        while abs(xvals[1] - xvals[2]) <= field_step / 2:
            k += 1
            xvals[2] = loop_fields[j - k] - B_offset
            yvals[2] = loop_moments[j - k] - M_offset

        ycal = y_polate(2, xvals, yvals, grid_fields[i])  #, grid_moments[i]
        grid_moments[i]=ycal
       # print(ycal)
    j -= 1
    grid_moments[n_loop // 2] = loop_moments[n_loop // 2] - M_offset

    for i in range(n_loop // 2 , n_loop):
        while not (loop_fields[j] - B_offset <= grid_fields[i]): # or (j <= n_loop // 2 + 1):
            if j <= n_loop // 2: #+ 1:
                break
            if j > n_loop // 2: #+ 1:
                j -= 1
        while not (loop_fields[j] - B_offset >= grid_fields[i]): # or (j >= n_loop):
            if j >= n_loop -1:
                break
            if j < n_loop -1:
                j += 1

        k = 1
        xvals[1] = loop_fields[j] - B_offset
        xvals[2] = loop_fields[j - k] - B_offset
        yvals[1] = loop_moments[j] - M_offset
        yvals[2] = loop_moments[j - k] - M_offset

        while abs(xvals[1] - xvals[2]) <= field_step / 2:
            k += 1
            xvals[2] = loop_fields[j - k] - B_offset
            yvals[2] = loop_moments[j - k] - M_offset

        ycal = y_polate(2, xvals, yvals, grid_fields[i])  #, grid_moments[i]
        grid_moments[i] = ycal
       # print(ycal)
    j = 2
    """   IRM loop fitting is fixed to linear for gridding purposes
    if polydegree > 1:
        for i in range(n2 + 1, n_loop // 2 - (n2 + 1) + 1):
            while not (loop_fields[j] - B_offset >= grid_fields[i]): # or (j <= 2):
                if j<=2:
                    break
                if j > 2:
                    j -= 1
            while not (loop_fields[j] - B_offset <= grid_fields[i]): # or (j >= n_loop // 2):
                if j >= n_loop // 2:
                    break
                if j < n_loop // 2:
                    j += 1

            for k in range(-n2, n2 + 1):
                x[k + n2 + 1] = loop_fields[j + k] - B_offset
                y[k + n2 + 1] = loop_moments[j + k] - M_offset

            # Call to polynomialfit1 function goes here
           # print(x)
            r2 = float()
            polynomialfit1(nsmooth, 0, polydegree, x, y) #polynomialfit1(nsmooth, 0, polydegree, x, y, y, r2)


            y1 = 0
            x1 = 1

            for k in range(polydegree + 1):
                if k > 0:
                    x1 *= grid_fields[i]
                y1 += y[k + 1] * x1

            grid_moments[i] = y1

        j -= 1

        for i in range(n_loop // 2 + n2 + 1, n_loop - (n2 + 1) + 1):
            while not (loop_fields[j] - B_offset <= grid_fields[i]) or (j <= n_loop // 2 + 1):
                if j > n_loop // 2 + 1:
                    j -= 1
            while not (loop_fields[j] - B_offset >= grid_fields[i]) or (j >= n_loop):
                if j < n_loop:
                    j += 1

            try:
                for k in range(-n2, n2 + 1):
                    x[k + n2 + 1] = loop_fields[j + k] - B_offset
                    y[k + n2 + 1] = loop_moments[j + k] - M_offset
            except:
                x1 = 0

            if j + k <= n_loop:
                # Call to polynomialfit1 function goes here
                # polynomialfit1(nsmooth, 0, polydegree, x, y, y, r2)
                # The implementation of polynomialfit1 function is needed

                y1 = 0
                x1 = 1

                for k in range(polydegree + 1):
                    if k > 0:
                        x1 *= grid_fields[i]
                    y1 += y[k + 1] * x1

                grid_moments[i] = y1
    """
    return grid_fields, grid_moments

def loop_test_linearity(n_looppoints, loop_fields, loop_moments):

    intercept, slope = linefit1(n_looppoints, loop_fields, loop_moments)
    n = n_looppoints // 2
    sst = 0
    ssr = 0
    SSD = 0
    SSLF = 0
    SSPE = 0
    ymean = 0
    ymean = sum(loop_moments) / (2 * n)

    for i in range(1, n+1):   # range(1, n + 1):
        b1 = loop_fields[i - 1]
        b2 = -loop_fields[i + n - 1]
        m1 = loop_moments[i - 1]
        m2 = -loop_moments[i + n - 1]
        mfit = b1 * slope + intercept
        sst += (m1 - ymean) ** 2
        ssr += (mfit - ymean) ** 2
        SSD += (mfit - m1) ** 2
        mfit = b2 * slope + intercept
        sst += (m2 - ymean) ** 2
        ssr += (mfit - ymean) ** 2
        SSD += (mfit - m2) ** 2
        # SSPE += ((m1 - m2) / 2) ** 2  Chat GPT botched the conversion
        SSPE += ((m1 - m2) ** 2) / 2

    SSLF = SSD - SSPE
    MSR = ssr
    MSD = SSD / (2 * n - 2)
    MSPE = SSPE / n
    MSLF = SSLF / (n - 2)
    FL = MSR / MSD
    FNL = MSLF / MSPE

    return FNL, slope, intercept

def loop_R2(n_looppoints, loop_fields, loop_moments, x26):
    n = n_looppoints
    min2 = 9E9
    max2 = -min2

    for i in range(n):
        loop_fields[i] -= x26

    for i in range(n):
        if loop_fields[i] > max2:
            max2 = loop_fields[i]

    for i in range(n):
        if loop_fields[i] < min2:
            min2 = loop_fields[i]

    min1 = -max2
    max1 = -min2
    n1 = n // 2
    i2 = 0
    #x1 = [None] * (n1 + 1)
    #y1 = [None] * (n1 + 1)
    x1 = np.ndarray([n1])
    y1 = np.ndarray([n1])
    x1[:]=0
    y1[:]=0


    n2 = 0
    for i in range(n1, n):
        x22 = loop_fields[i]   #changing x to x22
        if min1 < x22 < max1:
            while -loop_fields[i2] < x22:
                i2 += 1
            if i2 > 0:
                n2 += 1
                x1[n2] = loop_moments[i]
                dx = (-loop_fields[i2] - x22) / (-loop_fields[i2] + loop_fields[i2 - 1])
                dy = dx * (-loop_moments[i2] + loop_moments[i2 - 1])
                y = -loop_moments[i2] - dy
                y1[n2] = -y
    #V_shift,r2 = linefit1(n2, x1, y1)
    intercept, slope, rsqr = linefit2(n2,x1,y1)
   # print(intercept)
   # print(slope)
   # print(rsqr)
    V_shift = intercept / 2
    r2 = rsqr
    return r2, V_shift



def loop_R2v2(n_looppoints, loopfields, moments, x26):
    n = n_looppoints
    min2 = 9E9
    max2 = -min2
    loop_fields_Hshift = []
    for i in range(n):
        loop_fields_Hshift.append(loopfields[i] - x26)  #loop_fields.append(float(s1))

    for i in range(n):
        if loop_fields_Hshift[i] > max2:
            max2 = loop_fields_Hshift[i]

    for i in range(n):
        if loop_fields_Hshift[i] < min2:
            min2 = loop_fields_Hshift[i]

    min1 = -max2
    max1 = -min2
    n1 = n // 2
    i2 = 0
    x1 = np.zeros(n1)
    y1 = np.zeros(n1)

    n2 = 0
    for i in range(n1, n):
        x22 = loop_fields_Hshift[i]
        if min1 < x22 < max1:
            while -loop_fields_Hshift[i2] < x22:
                i2 += 1
            if i2 > 0:
                n2 += 1
                x1[n2] = moments[i]
                dx = (-loop_fields_Hshift[i2] - x22) / (-loop_fields_Hshift[i2] + loop_fields_Hshift[i2 - 1])
                dy = dx * (-moments[i2] + moments[i2 - 1])
                y = -moments[i2] - dy
                y1[n2] = -y

    intercept, slope, rsqr = linefit2(n2, x1, y1)
    V_shift = intercept / 2
    r2 = rsqr
    return r2, V_shift


def objective_function(x26):
    return loop_R2v2(n_looppoints, loop_fields, loop_moments, x26)[0]


def loop_Hshift_brent(n_looppoints, loopfields, moments, ax, bx, cx, tol):
   # def loop_R2(n_looppoints, loop_fields, loop_moments, H_shift, V_shift):
        # Define loop_R2 function here or import it from another module
   #     pass

    itmax = 100
    cgold = 0.3819660
    zeps = 1.0E-10

    a = ax if ax < cx else cx
    b = ax if ax > cx else cx
    v = bx
    w = v
    x26 = v   #changing x to x26 because seems like python does not like me using x
    e = 0.0
    #fx = -loop_R2(n_looppoints, loop_fields, loop_moments, x)
    fx=loop_R2v2(n_looppoints, loopfields, moments, x26)[0]
    fx = fx * -1
    #print(fx)
    fv = fx
    fw = fx
    d = 0
    goto_1 = True
    goto_2 = True
    for iter in range(1, itmax + 1):
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x26) + zeps
        tol2 = 2.0 * tol1

        if abs(x26 - xm) <= (tol2 - 0.5 * (b - a)):
            break

        if abs(e) > tol1:
            r = (x26 - w) * (fx - fv)
            q = (x26 - v) * (fx - fw)
            p = (x26 - v) * q - (x26 - w) * r
            q = 2.0 * (q - r)

            if q > 0.0:
                p = -p
            q = abs(q)
            etemp = e
            e = d
            goto_1 = True

            if abs(p) >= abs(0.5 * q * etemp) or p <= q * (a - x26) or p >= q * (b - x26):
                goto_1 = True
            else:
                d = p / q
                u = x26 + d

                if (u - a) < tol2 or (b - u) < tol2:
                    d = sign(tol1, xm - x26)
                goto_2 = True
                goto_1 = False

        if goto_1: #or (not goto_1 and abs(d) >= tol1):
            if x26 >= xm:
                e = a - x26
            else:
                e = b - x26
            d = cgold * e

        if goto_2: #or (not goto_2 and abs(d) >= tol1):
            if abs(d) >= tol1:
                u = x26 + d
            else:
                u = x26 + sign(tol1, d)
            fu=loop_R2v2(n_looppoints, loopfields, moments, u)[0]
            fu = fu * -1
            if fu <= fx:
                if u >= x26:
                    a = x26
                else:
                    b = x26
                v = w
                fv = fw
                w = x26
                fw = fx
                x26 = u
                fx = fu
            else:
                if u < x26:
                    a = u
                else:
                    b = u
                if fu <= fw or w == x26:
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                elif fu <= fv or v == x26 or v == 2:
                    v = u
                    fv = fu

    boff1 = x26
    r2 = fx
    #moff1 = loop_R2v2(n_looppoints, loopfields, moments, u)[1]
    moff1 = loop_R2v2(n_looppoints, loop_fields, loop_moments, u)[1]
    return r2, boff1, moff1

def sign(a, b):
    return abs(a) if b > 0.0 else -abs(a)

#import math

def loop_errorcal(n_looppoints, loopfields, moments):
   # def loop_Hshift_brent(n_looppoints, loopfields, moments, ax, bx, cx, tol, xmin, y0):
        # Implement loop_Hshift_brent here or import it from another module
   #     pass

    #x = [None] * n_looppoints
    #y = [None] * n_looppoints
    #r2, need to get r2 somehow
    r2, boff1, moff1 = loop_Hshift_brent(n_looppoints, loopfields, moments, -loopfields[0] / 2, 0, loopfields[0] / 2, 1E-6)
    #result = brent(objective_function, brack=(-loopfields[0] / 2, 0, loopfields[0] / 2), tol=1E-6)


    M_offset = moff1
    B_offset = boff1

    if 1 + r2 > 0:
        M_sn = math.sqrt(1 + r2)  # =sqrt(1-R^2) = sqrt(SSD/SST) ~ noise/signal
    else:
        M_sn = 0

    return M_offset, B_offset, M_sn


def loop_errorcal2(n_looppoints, loop_fields, loop_moments, M_offset, B_offset):
    n = n_looppoints
    min2 = 9E9
    max2 = -min2

    for i in range(n):
        loop_fields[i] -= B_offset

    for i in range(n):
        if loop_fields[i] > max2:
            max2 = loop_fields[i]

    for i in range(n):
        if loop_fields[i] < min2:
            min2 = loop_fields[i]

    min1 = -max2
    max1 = -min2
    n1 = n // 2
    i2 = 0
    n2 = 0
    ErrX = []
    ErrY = []
    for i in range(n1, n):
        x = loop_fields[i]
        y0 = loop_moments[i]

        if (x - min1 > 1E-10) and (x < max1):
            while -loop_fields[i2] < x:
                i2 += 1

            if i2 > 0:
                #n2 += 1   if this is to be a 0 bassed array, need to move inc to end or use append
                ErrX.append(-x)    #ErrX[n2] = -x
                dx = (-loop_fields[i2] - x) / (-loop_fields[i2] + loop_fields[i2 - 1])
                dy = dx * (-loop_moments[i2] + loop_moments[i2 - 1])
                y = -loop_moments[i2] - dy
                ErrY.append(y0 - y)       #ErrY[n2] = y0 - y


    return ErrX, ErrY


def loop_delta_M(n_looppoints, loopfields, moments):
    # Initialize variables
    n2 = n_looppoints // 2
    xarr = np.zeros(n2 + 1)  # array needs to be filled with zeros
    yarr = np.zeros(n2 + 1)
    mrh = np.zeros(n2 + 1)

    E_hys = 0
    Mrhmax = 0
    noise1 = 0

    # Main loop
    for i in range(1, n2 + 1):
        j = n_looppoints + 51 - i
        if j > n_looppoints:
            j = n_looppoints
        while not (loopfields[j - 1] <= loopfields[i - 1] or j == n_looppoints // 2):
            j -= 1

        xarr[0] = loopfields[j - 1]
        xarr[1] = loopfields[j]
        yarr[0] = moments[j - 1]
        yarr[1] = moments[j]

        if xarr[0] == xarr[1]:
            y1 = yarr[0]
        else:
            y1 = y_polate(2, xarr, yarr, loopfields[i - 1])

        mrh[i - 1] = (moments[i - 1] - y1) / 2

        if mrh[i - 1] > Mrhmax:
            Mrhmax = mrh[i - 1]

        if i > 1:
            E_hys -= (mrh[i - 1] + mrh[i - 2]) / 2 * (loopfields[i - 1] - loopfields[i - 2])
            noise1 += (mrh[i - 1] - mrh[i - 2]) ** 2

    noise1 = np.sqrt(noise1 / n2)
    drift1 = np.sum(mrh[:10])

    # Calculate Brh, Bih
    MS = np.mean(moments[2:6])
    j = 0
    while not (moments[j] - mrh[j] <= MS / 2 or j == n_looppoints // 2):
        j += 1
    xarr[0] = loopfields[j - 1]
    xarr[1] = loopfields[j - 2]
    yarr[0] = moments[j - 1] - mrh[j - 1]
    yarr[1] = moments[j - 2] - mrh[j - 2]
    x1 = x_polate(2, xarr, yarr, MS / 2)
    Bih = x1

    j = 0
    while not (mrh[j] >= Mrhmax / 2 or j == n_looppoints // 2):
        j += 1
    xarr[0] = loopfields[j - 1]
    xarr[1] = loopfields[j - 2]
    yarr[0] = mrh[j - 1]
    yarr[1] = mrh[j - 2]
    x1 = x_polate(2, xarr, yarr, Mrhmax / 2)
    Brh = x1

    return mrh, E_hys, Brh, Bih

def loop_drift(n_looppoints, loop_fields, loop_moments, err_moments):
    n2 = n_looppoints // 2
    n4 = n_looppoints // 4
    x1 = np.zeros(n2 + 1)
    y1 = np.array(err_moments[:n2 + 1])

    n = 5  # smooth 2n+1 running poly
    n1 = 2  # poly deg
    x2 = np.zeros(2 * n + 2)
    y2 = np.zeros(2 * n + 2)

    for i in range(n + 1, n2 - n -2):  #for i in range(n + 1, n2 - n + 1):
        for j in range(-n, n + 1):
            jplusn = j + n
            iplusj = i + j
            x2[j + n] = loop_fields[i + j]
            y2[j + n] = err_moments[i + j]
        y0 = polynomialfit2(2 * n + 1, 0, n1, x2, y2, x2[n])
        y1[i] = y0  # y1 is smoothed error curve

    s1 = np.sum(np.square(y1)) / n2

    s2 = 0
    r2 = 0
    for i in range(n2-2):
        if abs(y1[i]) > s2:
            s2 = abs(y1[i])
            r2 = loop_fields[i]

    grid_moments_cor = np.copy(loop_moments)
    corr_applied = abs(r2) > loop_fields[0] * 0.75
    if corr_applied:
        for i in range(1, n4 + 1):
            grid_moments_cor[i - 1] -= y1[i]
            grid_moments_cor[n_looppoints - i] -= y1[n2 - i]

    return grid_moments_cor, corr_applied

def loop_driftcorr2(n_looppoints, loop_fields, loop_moments, err_moments):
    n = n_looppoints // 2

    for i in range(0, n_looppoints):   #for i in range(1, n_looppoints + 1):
        grid_moments_cor[i] = loop_moments[i]
    k = 3
    for i in range(0 + k, n - k -2):    #for i in range(1 + k, n - k + 1):
        y = 0
        for j in range(-k, k + 1):
            y += err_moments[i + j]
        y /= (2 * k + 1)
        grid_moments_cor[i] = loop_moments[i] - y
    corr_applied = True

    return grid_moments_cor, corr_applied

def loop_test_saturation(n_looppoints, loop_fields, loop_moments, err_moments):
    """
    def betai(a, b, x):
        # Define betai function here or import from another module
        pass

    def linefit1(n, x, y, intercept, slope, r2):
        # Define linefit1 function here or import from another module
        pass
    """
    FNS60 = 0.0
    FNS70 = 0.0
    FNS80 = 0.0
    MS = 0.0
    Xhf = 0.0

    x1 = 0.97
    x0 = 0.0
    xcutoff2 = x1 * loop_fields[0]
    n_hfpar = 4
    n = 0
    bmax = 0.0

    #xarr = np.zeros(n_looppoints)
    #yarr = np.zeros(n_looppoints)


    for j in range(1, 4):
        xarr = []
        yarr = []
        if j == 1:
            xcutoff1 = 0.6 * loop_fields[0]
        elif j == 2:
            xcutoff1 = 0.7 * loop_fields[0]
        elif j == 3:
            xcutoff1 = 0.8 * loop_fields[0]

        n = 0   #n = 0
        bmax = 0.0

        for i in range(n_looppoints):
            if abs(loop_fields[i]) >= xcutoff1 and abs(loop_fields[i]) <= xcutoff2:
                #n += 1
                xarr.append(abs(loop_fields[i]))
                yarr.append(np.sign(loop_fields[i]) * loop_moments[i])
                n += 1

            if abs(loop_fields[i]) > bmax:
                bmax = abs(loop_fields[i])

        #for idx, val in enumerate(yarr):
            #print(f"Index {idx}: {val}")


        intercept, slope, r2 = linefit2(n, xarr, yarr)   #intercept, slope, r2 = linefit2(n, np.array(x), np.array(y))

        k = n // 2
        sst = 0.0
        ssr = 0.0
        SSD = 0.0
        SSLF = 0.0
        SSPE = 0.0
        ymean = np.mean(yarr)
        """
        for i in range(1, n+1):
            sst += (y[i-1] - ymean) ** 2
        """
        for i in range(1, k+1):
            b1 = xarr[i-1]
            b2 = xarr[i + k - 1]
            m1 = yarr[i-1]
            m2 = yarr[i + k - 1]

            mfit = b1 * slope + intercept
            sst += (m1 - ymean) ** 2
            ssr += (mfit - ymean) ** 2
            SSD += (mfit - m1) ** 2

            mfit = b2 * slope + intercept
            sst += (m2 - ymean) ** 2
            ssr += (mfit - ymean) ** 2
            SSD += (mfit - m2) ** 2
            SSPE += (m1 - m2)**2/ 2

        SSLF = SSD - SSPE
        MSR = ssr
        MSD = SSD / (n - 2)
        MSPE = SSPE / k

        if (n - 2 - k) > 0:
            MSLF = SSLF / (n - 2 - k)
        else:
            MSLF = 0.0

        FL = MSR / MSD

        if j == 1:
            FNS60 = MSLF / MSPE
            df1 = n - 2 - k
            df2 = k
            p = df1 / (df1 + df2 / FNS60)
            if 0 <= p <= 1:
                p60 = sc.betainc(df1 / 2, df2 / 2, p)  #Regularized incomplete beta function.
            else:
                p60 = 0.0
        elif j == 2:
            FNS70 = MSLF / MSPE
            df1 = n - 2 - k
            df2 = k
            p = df1 / (df1 + df2 / FNS70)
            if 0 <= p <= 1:
                p70 = sc.betainc(df1 / 2, df2 / 2, p)  #Regularized incomplete beta function.
            else:
                p70 = 0.0
        elif j == 3:
            FNS80 = MSLF / MSPE
            df1 = n - 2 - k
            df2 = k
            if FNS80 > 0:
                p = df1 / (df1 + df2 / FNS80)
                if 0 <= p <= 1:
                    p80 = sc.betainc(df1 / 2, df2 / 2, p)   #Regularized incomplete beta function.
                else:
                    p80 = 0.0
            else:
                p80 = 0.0

    if FNS80 < 2.5:  #saturated at 80%
        x0 = 0.8
    if FNS70 < 2.5:  #saturated at 70%
        x0 = 0.7
    if FNS60 < 2.5:  #saturated at 60%
        x0 = 0.6
    if x0 < 0.5:   #nonsaturated linear fit
        x0 = 0.92
        x1 = 0.97
        MS = 0
        Xhf = 0
        r2 = 0
    else:
        j = 0
        i = n_looppoints // 2
        xarr = np.zeros(n_looppoints // 4)
        yarr = np.zeros(n_looppoints // 4)
        while abs(loop_fields[i]) <= x1 * loop_fields[0]:
            i += 1
            if abs(loop_fields[i]) <= x1 * loop_fields[0]:
                j += 1
                xarr[j-1] = loop_fields[i-1]
                yarr[j-1] = loop_moments[i-1]

        MS, Xhf = linefit1(j, xarr, yarr)     #MS, Xhf, r2 = linefit2(j, x, y)
        MS = -MS
    return FNS60, FNS70, FNS80, MS, Xhf

def loop_fit_saturation(n_looppoints, loop_fields, loop_moments, fixed_beta):
    # nonlinear fit
    npoints = 0
    x0 = 0.6
    x1 = 0.95
    bmax = loop_fields[0]
    n_hfpar = 3
    MHData = []
    dmx = []

    for j in range(n_looppoints // 4):
        if abs(loop_fields[j]) / bmax > x0 and abs(loop_fields[j]) / bmax < x1:
            npoints += 1
            MHData.append([loop_fields[j], loop_moments[j]] if loop_fields[j] > 0 else [-loop_fields[j], -loop_moments[j]])
            dmx.append([loop_fields[j], 1])
            if n_hfpar > 2:
                dmx[-1].append(1 / dmx[-1][0])  #dmx[-1].append(1 / dmx[-1][1])
            if n_hfpar > 3:
                dmx[-1].append(1 / (dmx[-1][0] * dmx[-1][0]))

    xarr = np.array([data[0] for data in MHData])
    yarr = np.array([data[1] for data in MHData])

    #intercept, slope, _ = np.linalg.lstsq(np.vstack([x, np.ones(len(x))]).T, y, rcond=None)
    #intercept, slope, _ = np.linalg.lstsq(np.vstack([x, np.ones(len(x))]).T, y)
    intercept, slope, r2lin = linefit2(npoints, xarr, yarr);

    # normalize by measurement std error
    #M_err = 1.0e-7  # adjust as needed
    for i in range(len(MHData)):
        for i1 in range(n_hfpar):
            dmx[i][i1] /= M_err
        MHData[i][1] /= M_err
   # SV_U, SV_A, SV_V = svd(np.array(dmx), full_matrices=False)
   # SV_cov = np.linalg.inv(SV_V.T @ np.diag(svdvals(SV_A)) ** 2 @ SV_V)

    beta_ = np.zeros(35)
    beta_Ms = np.zeros(35)
    beta_Xhf = np.zeros(35)
    rsquare = np.zeros(35)
    #MHcoeffs = np.zeros((5, 1))
    dmxt = np.zeros((n_hfpar, npoints))
    SV_A = np.zeros((npoints, n_hfpar))

    if fixed_beta:
        nbeta = 1
    else:
        nbeta = 11

    z3 = -1E10 #Fabian fit??
    m_measured = np.zeros(npoints)
    m_fit = np.zeros(npoints)


    for i in range(nbeta):
        dmxtdmx1 = np.zeros((3,3))
        if fixed_beta:
            beta1 = beta_fixed
        else:
            beta1 = -(i + nbeta - 2) / (nbeta - 1)  # -1.0 to -2.0

        for j in range(npoints):
            dmx[j][0] = 1
            dmx[j][1] = MHData[j][0]
            dmx[j][2] = pow(MHData[j][0], beta1)

         # transpose dmx
        for i1 in range(0, n_hfpar):
            for j in range(0, npoints):
                dmxt[i1][j] = dmx[j][i1]
            for j in range(0,npoints):
                SV_A[j][i1] = dmx[j][i1]

        """
        dmx_temp = [[1, data[0], np.exp(data[0] * beta1)] for data in MHData]
        for j in range(len(dmx_temp)):
            for k in range(0,n_hfpar):
                dmx_temp[j][k] /= M_err

            #dmx_temp[j][0] /= M_err
            #dmx_temp[j][1] /= M_err
            #dmx_temp[j][2] /= M_err
        for i in range(len(dmx_temp)):
            print(dmx_temp[i])
        sv_b = np.array([data[1] for data in MHData])
        """

        """
        _, sv_w, SV_V = svd(np.array(dmx_temp), full_matrices=False)
        sv_x = solve(np.diag(svdvals(SV_A))[:, None] * SV_V, sv_b)
        SV_cov = np.linalg.inv(SV_V.T @ np.diag(svdvals(SV_A)) ** 2 @ SV_V)
        """

        dmxtdmx = np.matmul(dmxt, dmx);

        for i1 in range(3):
            for j in range(3):
                dmxtdmx1[i1][j] = dmxtdmx[i1][j]

        dmxtdmxinv = np.linalg.inv(dmxtdmx1)
        dmxtdmx1 = np.matmul(dmxt, MHData)
        MHcoeffs = np.matmul(dmxtdmxinv, dmxtdmx1)
        fitdata = np.matmul(dmx, MHcoeffs)

        """
        MHcoeffs = solve(np.array(dmxt).T @ np.array(dmxt), np.array(dmxt).T @ sv_b)
        fitdata = np.array(dmx_temp) @ MHcoeffs
        """
        z = 0
        z0 = 0
        z1 = 0
        z2 = 0
        for j in range(npoints):
            z0 += MHData[j][1]
        z0 /= npoints
        for j in range(npoints):
            z += (MHData[j][1]- z0)**2
            z1 += (fitdata[j][1]-z0)**2
            z2 += (fitdata[j][1]-MHData[j][1])**2
            m_measured[j] = MHData[j][1]
            m_fit[j] = fitdata[j][1]

        x1, y1, r22 = linefit2(npoints, m_measured, m_fit)


        """
        z = np.sum((sv_b - np.mean(sv_b)) ** 2)
        z1 = np.sum((fitdata - np.mean(sv_b)) ** 2)
        z2 = np.sum((fitdata - sv_b) ** 2)
        m_measured = sv_b
        m_fit = fitdata
        """

        beta_[i] = beta1
        beta_Ms[i] = MHcoeffs[1][1] * M_err
        beta_Xhf[i] = MHcoeffs[2][1] * M_err
        rsquare[i] = r22

        if rsquare[i] > z3:  # Best fit
            MS = MHcoeffs[0][1] * M_err
            Xhf = MHcoeffs[1][1] * M_err
            alpha = MHcoeffs[2][1] * M_err
            beta = beta1
            z3 = rsquare[i]
            F_lnl = (z3 - r2lin) / (1 - z3) * (npoints - 4) / 2

    if (alpha > 0 or F_lnl < 3) and not fixed_beta:
        SV_A = np.array(dmx)[:, :2]
        sv_b = np.array([data[1] for data in MHData])
        _, sv_w, SV_V = svd(SV_A, full_matrices=False)
        sv_x = solve(np.diag(svdvals(SV_A))[:, None] * SV_V, sv_b)
        SV_cov = np.linalg.inv(SV_V.T @ np.diag(svdvals(SV_A)) ** 2 @ SV_V)
        MS = sv_x[0] * M_err
        Xhf = sv_x[1] * M_err
        alpha = 0

    return  MS, Xhf, alpha, beta, F_lnl

def loop_errorcal2(n_looppoints, loop_fields, loop_moments, M_offset, B_offset):
    n = n_looppoints
    min2 = 9E9
    max2 = -min2
    ErrX = np.zeros(n_looppoints)
    ErrY = np.zeros(n_looppoints)
    for i in range(n):
        loop_fields[i] -= B_offset

    for i in range(n):
        if loop_fields[i] > max2:
            max2 = loop_fields[i]

    for i in range(n):
        if loop_fields[i] < min2:
            min2 = loop_fields[i]

    min1 = -max2
    max1 = -min2
    n1 = n // 2
    i2 = 1
    n2 = 0

    for i in range(n1, n):
        x = loop_fields[i]
        y0 = loop_moments[i]
        if (x - min1 > 1E-10) and (x < max1):
            while -loop_fields[i2] < x:
                i2 += 1
            if i2 > 0:
                n2 += 1
                ErrX[n2] = -x
                dx = (-loop_fields[i2] - x) / (-loop_fields[i2] + loop_fields[i2 - 1])
                dy = dx * (-loop_moments[i2] + loop_moments[i2 - 1])
                y = -loop_moments[i2] - dy
                ErrY[n2] = y0 - y
    return ErrX, ErrY


def loop_intercepts(n_looppoints, loop_fields, loop_moments):
    xarr = np.zeros(8)
    yarr = np.zeros(8)
    aarr = np.zeros(8)
    Hcminus = 0
    Hcplus = 0
    Mrplus = 0
    Mrminus = 0
    XMr = 0
    XHc = 0

    i = 0
    while loop_fields[i] > 0:
        i += 1

    for j in range(0, 8):
        xarr[j] = loop_fields[i + 4 - j]
        yarr[j] = loop_moments[i + 4 - j]

    # quadratic fit near H=0
    a = np.polyfit(xarr, yarr, 2)
    Mrplus = a[2]
    XMr = a[1]

    while loop_moments[i] > 0:
        i += 1

    for j in range(0, 8):
        xarr[j] = loop_fields[i + 4 - j] - loop_fields[i]
        yarr[j] = loop_moments[i + 4 - j]

    # cubic fit near M=0
    a = np.polyfit(xarr, yarr, 3)
    Hcminus = loop_fields[i] - a[3] / a[2]
    XHc = a[2]

    i = n_looppoints // 2
    while loop_fields[i] < 0:
        i += 1

    for j in range(0, 8):
        xarr[j] = loop_fields[i + 4 - j]
        yarr[j] = loop_moments[i + 4 - j]

    a = np.polyfit(xarr, yarr, 2)
    Mrminus = a[2]
    XMr = (XMr + a[1]) / 2

    while loop_moments[i] < 0 and i < n_looppoints - 3:
        i += 1

    for j in range(0, 8):
        xarr[j] = loop_fields[i + 4 - j] - loop_fields[i]
        yarr[j] = loop_moments[i + 4 - j]

    a = np.polyfit(xarr, yarr, 3)
    Hcplus = loop_fields[i] - a[3] / a[2]
    XHc = (XHc + a[2]) / 2

    return Hcminus, Hcplus, Mrplus, Mrminus, XMr, XHc

def simple_HF_fit(n_looppoints, loop_fields, loop_moments, pctsat):
    x1 = 0.97
    x0 = pctsat
    j = 0
    i = n_looppoints // 2
    xarr = np.zeros(n_looppoints // 4)
    yarr = np.zeros(n_looppoints // 4)
    while abs(loop_fields[i]) >= x0 * loop_fields[0]:
        i += 1
        if abs(loop_fields[i]) <= x1 * loop_fields[0]:
            j += 1
            xarr[j-1] = loop_fields[i-1]
            yarr[j-1] = loop_moments[i-1]
    MS, Xhf = linefit1(j, xarr, yarr)     #MS, Xhf, r2 = linefit2(j, x, y)
    MS = -MS

    return MS, Xhf

#  main  -------------------------------------------------------------------
mu0 = 1.2566370614E-6
loop_fields = []
loop_moments = []
f=open(r"G:\My Drive\python\ned15c_field.txt")
for x in f:
    s1=x
    if len(s1) >2:loop_fields.append(float(s1))
f.close()
f=open(r"G:\My Drive\python\ned15c_moment.txt")
for x in f:
    s1=x
    if len(s1) >2:loop_moments.append(float(s1)/.000153)
f.close()
n_loop = len(loop_fields)
polydegree = 1
nsmooth = 3

B_offset = 0
M_offset = 0
acheckvar= 0

#print(loop_moments)

# first call to grid routine ignores offsets for linearity test
grid_fields, grid_moments = loop_grid(n_loop, polydegree, nsmooth, loop_fields, loop_moments, B_offset, M_offset)
print("grid moments[0] = ",grid_moments[0])
print("grid fields[0] = ", grid_fields[0])

moment_sub = branch_sub(n_loop, grid_moments)

#plt.subplot(1,3,1)
#plt.plot(loop_fields, loop_moments)

#plt.subplot(1,3,2)
#plt.plot(grid_fields, grid_moments)




#xpoints = grid_fields[n_loop // 2:]

#plt.subplot(1,3,3)
#ypoints = moment_sub
#plt.plot(xpoints, ypoints)
#plt.show()

n_looppoints = n_loop

#Lets start with a simple linerar fit to the high field data above 80% of the maximum field
pctsat = 0.8
MS, Xhf = simple_HF_fit(n_looppoints, loop_fields, loop_moments, pctsat)
print("MS = ", MS)
print("Xhf = ", Xhf)


#subrtact Xhf from loop
loop_moments_f = np.zeros(n_looppoints)
for i in range(n_looppoints):
    loop_moments_f[i] = grid_moments[i] - Xhf * grid_fields[i]   #slope correction

plt.subplot(1,2,1)
plt.plot(grid_fields, grid_moments)
plt.subplot(1,2,2)
plt.plot(grid_fields, loop_moments_f)
plt.show()


FNL, slope, intercept = loop_test_linearity(n_looppoints, grid_fields, grid_moments)
print("FNL = ",FNL)
print("slope = ",slope)
print("intercept = ",intercept)

M_offset, B_offset, M_sn = loop_errorcal(n_looppoints, grid_fields, grid_moments);
print('good luck')
print("M_offset = ",M_offset)
print("B_offset = ",B_offset)

"""
A quality factor Q is calculated for both the whole curve, Q, and the ferromagnetic part, Qf.
Q is the decimal log of the signal/noise ratio, calculated from the mean square mismatch
between symmetrically equivalent moments, and typically ranges between 0 (s/n~1) and ~3 (s/n~1000).
A higher value indicates better quality data.

"""


if M_sn > 0:
    Q = math.log(M_sn, 10) * -1
else:
    Q = 10
print("Q = ",Q)

# second call to grid routine removes offsets
grid_fields, grid_moments = loop_grid(n_loop, polydegree, nsmooth, loop_fields, loop_moments, B_offset, M_offset)


#plt.subplot(1,3,1)
#plt.plot(loop_fields, loop_moments)

#plt.subplot(1,3,2)
#plt.plot(grid_fields, grid_moments)

#xpoints = grid_fields[n_loop // 2:]

#plt.subplot(1,3,3)
#ypoints = moment_sub
#plt.plot(xpoints, ypoints)
#plt.show()

print ("grid_moments = ",grid_moments[0])
print ("grid_fields = ",grid_fields[0])

ErrX, ErrY = loop_errorcal2(n_looppoints, grid_fields, grid_moments, 0, 0)
for i in range(10):
    print("error x and y = ",ErrX[i],", ",ErrY[i])


M_err = 0
for i in range(0, n_looppoints//2 -2):
    M_err += (ErrY[i])**2
M_err = math.sqrt(M_err / (n_looppoints//2 ) )
print("n_looppoints = ", n_looppoints//2)
print("M_err = ", M_err)

mrh, E_hys, Brh, Bih = loop_delta_M(n_looppoints, grid_fields, grid_moments)
print('loop_delta_m first go')
#print("mrh =", mrh)
print("E_hys =", E_hys)
print("Brh =", Brh)
print("Bih =", Bih)

"""
plt.subplot(1,4,1)
plt.plot(loop_fields, loop_moments)

plt.subplot(1,4,2)
plt.plot(grid_fields, grid_moments)

xpoints = grid_fields[n_loop // 2:]
apoints = grid_fields[n_loop // 2:]
plt.subplot(1,4,3)
ypoints = moment_sub
plt.plot(xpoints, ypoints)

plt.subplot(1,4,4)
zpoints = mrh[:-1]
plt.plot(apoints, zpoints)
plt.show()
"""
grid_moments_cor, corr_applied = loop_drift(n_looppoints, grid_fields, grid_moments, ErrY)
#print(grid_moments_cor[10])
print("corr_applied = ", corr_applied)

if not corr_applied:
    grid_moments_cor, corr_applied = loop_driftcorr2(n_looppoints, grid_fields, grid_moments, ErrY)

print(grid_moments_cor[10])
print("corr_applied = ",corr_applied)
"""
plt.subplot(1,2,1)
plt.plot(grid_fields, grid_moments)
plt.subplot(1,2,2)
plt.plot(grid_fields, grid_moments_cor)
plt.show()
"""
if corr_applied:
    mrh, E_hys, Brh, Bih = loop_delta_M(n_looppoints, grid_fields, grid_moments_cor)
print('loop_delta_m second go')
#print('mrh = ',mrh)
print('E_hys = ',E_hys)
print('Brh = ',Brh)
print('Bih = ',Bih)

FNS60, FNS70, FNS80, MS, Xhf = loop_test_saturation(n_looppoints, grid_fields, grid_moments_cor, ErrY)

print('FNS80 = ',FNS80)
print('FNS70 = ',FNS70)
print('FNS60 = ',FNS60)
print('MS = ',MS)
print('Xhf = ',Xhf)

if FNS70 > 2.5 and FNS80 > 2.5 and FNS60 > 2.5:
    fixed_beta = False  #use fixed beta for anisotropy or otherwise let beta float for indivitual loops
    Ms, Xhf, alphaHF, betaHF, F_nlfit = loop_fit_saturation(n_looppoints, grid_fields, grid_moments_cor, fixed_beta)
    print('MS = ',Ms)
    print('Xhf = ',Xhf*mu0)
    print('alphaHF = ',alphaHF)
    print('betaHF = ',betaHF)
    print('F_nlfit = ',F_nlfit)
else:
    if FNS60 < 2.5:
        percentage = 60
    elif FNS70 < 2.5:
        percentage = 70
    else:
        percentage = 80

loop_moments_f = np.zeros(n_looppoints)
for i in range(n_looppoints):
    loop_moments_f[i] = grid_moments_cor[i] - Xhf * grid_fields[i]   #slope correction

plt.subplot(1,2,1)
plt.plot(grid_fields, grid_moments_cor)
plt.subplot(1,2,2)
plt.plot(grid_fields, loop_moments_f)
plt.show()




moffset = 0
Boffset = 0
x1arr = np.zeros(n_looppoints)
y1arr = np.zeros(n_looppoints)
ErrX, ErrY = loop_errorcal2(n_looppoints, grid_fields, loop_moments_f, moffset,Boffset)

x1 = 0
y1 = 0
for i in range(n_looppoints//2):
    x1 += ErrY[i]**2
    y1 += loop_moments_f[i]**2

M_sn = math.sqrt(x1 / y1)
if M_sn > 0:
    Qf = math.log(M_sn, 10) * -1
else:
    Qf = 10
print("Qf = ", Qf)





Hcminus, Hcplus, Mrplus, Mrminus, XMr, XHc = loop_intercepts(n_looppoints, grid_fields, grid_moments)

print('Hcminus = ',Hcminus)
print('Hcplus = ',Hcplus)
print('Mrplus = ',Mrplus)
print('Mrminus = ',Mrminus)
print('XMr = ',XMr)
print('XHc = ',XHc)
print('Xmr = ',XMr*mu0)
print('-------------------------------------')

Hcminus, Hcplus, Mrplus, Mrminus, X1, X2 = loop_intercepts(n_looppoints, grid_fields, loop_moments_f)
print('Hcminus = ',Hcminus)
print('Hcplus = ',Hcplus)
print('Mrplus = ',Mrplus)
print('Mrminus = ',Mrminus)
print('XMr = ',X1)
print('XHc = ',X2)

sigma_hys = E_hys / (2 * Hcplus * Ms)

if sigma_hys > 0:
    sigma_hys = math.log(sigma_hys, 10) * -1
    print('sigma_hys = ',sigma_hys)
else:
    print('sigma_hys = N/A')
    sigma_hys = -90
