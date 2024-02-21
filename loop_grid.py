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
from numpy.linalg import svd, lstsq

def y_polate(n, xarr, yarr, xval, yval):
    # linear interpolation/extrapolation, returns best-fit y for a specified x (e.g. for Mr)
    intercept, slope = linefit1(n, xarr, yarr)
    ycal = intercept + slope * xval

    return ycal

def linefit1(n, xarr, yarr):
    sum_x = sum(xarr)
    sum_y = sum(yarr)
    sum_xy = sum(x * y for x, y in zip(xarr, yarr))
    sum_x_squared = sum(x ** 2 for x in xarr)

    slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x ** 2)
    intercept = (sum_y - slope * sum_x) / n

    return intercept, slope

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

def branch_sub(n_loop, grid_moments):

    j = n_loop // 2
    moment_sub = np.ndarray([j])
    for i in range (j):
        moment_sub[i] = grid_moments[i]-grid_moments[n_loop-1-i]
    return moment_sub


def loop_grid(n_loop, polydegree, nsmooth, loop_fields, loop_moments, B_offset, M_offset):

    n2 = nsmooth // 2  # make nsmooth odd; n2=# on each side, + 1 in center
    xyada = np.ndarray([nsmooth+3])
    yyada = np.ndarray([nsmooth+3])
    #x = [0.0] * (nsmooth + 3)
    #y = [0.0] * (nsmooth + 3)

    grid_fields = np.ndarray([n_loop])
    grid_moments = np.ndarray([n_loop])
    #grid_fields = [0.0] * (n_loop + 1)
    #grid_moments = [0.0] * (n_loop + 1)



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
        while not (loop_fields[j] - B_offset <= grid_fields[i]) or (j >= n_loop // 2):
            if j < n_loop // 2:
                j += 1

        k = 1
        xyada[1] = loop_fields[j] - B_offset
        xyada[2] = loop_fields[j - k] - B_offset
        yyada[1] = loop_moments[j] - M_offset
        yyada[2] = loop_moments[j - k] - M_offset

        while abs(xyada[1] - xyada[2]) <= field_step / 2:
            k += 1
            xyada[2] = loop_fields[j - k] - B_offset
            yyada[2] = loop_moments[j - k] - M_offset

        ycal = y_polate(2, xyada, yyada, grid_fields[i], grid_moments[i])
        grid_moments[i]=ycal
        #print(grid_moments[i])
    j -= 1
    grid_moments[n_loop // 2] = loop_moments[n_loop // 2] - M_offset

    for i in range(n_loop // 2, n_loop):
        while not (loop_fields[j] - B_offset <= grid_fields[i]): # or (j <= n_loop // 2 + 1):
            if j <= n_loop // 2 + 1:
                break
            if j > n_loop // 2 + 1:
                j -= 1
        while not (loop_fields[j] - B_offset >= grid_fields[i]) or (j >= n_loop):
            if j < n_loop:
                j += 1

        k = 1
        xyada[1] = loop_fields[j] - B_offset
        xyada[2] = loop_fields[j - k] - B_offset
        yyada[1] = loop_moments[j] - M_offset
        yyada[2] = loop_moments[j - k] - M_offset

        while abs(xyada[1] - xyada[2]) <= field_step / 2:
            k += 1
            xyada[2] = loop_fields[j - k] - B_offset
            yyada[2] = loop_moments[j - k] - M_offset

        ycal = y_polate(2, xyada, yyada, grid_fields[i], grid_moments[i])
        grid_moments[i] = ycal
        #print(grid_moments[i])
    j = 2

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
            print(x)
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

    return grid_fields, grid_moments

# Example usage

loop_fields = []
loop_moments = []
f=open("c:\glop3_field.txt")
for x in f:
    s1=x
    if len(s1) >2:loop_fields.append(float(s1))
f.close()
f=open("c:\glop3_moment.txt")
for x in f:
    s1=x
    if len(s1) >2:loop_moments.append(float(s1))
f.close()
n_loop = 1442
polydegree = 1
nsmooth = 3
#loop_fields = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
#loop_moments = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
B_offset = 0
M_offset = 0
acheckvar= 0
#grid_fields = [0.0] * (n_loop + 1)
#grid_moments = [0.0] * (n_loop + 1)


grid_fields, grid_moments = loop_grid(n_loop, polydegree, nsmooth, loop_fields, loop_moments, B_offset, M_offset)

moment_sub = branch_sub(n_loop, grid_moments)

plt.subplot(1,3,1)
plt.plot(loop_fields, loop_moments)

plt.subplot(1,3,2)
plt.plot(grid_fields, grid_moments)

xpoints = grid_fields[721:]

plt.subplot(1,3,3)
ypoints = moment_sub
plt.plot(xpoints, ypoints)
plt.show()
