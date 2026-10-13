from sympy import symbols, diff, sqrt

x, y, px, py, xt, yt, xm, ym, xl, yl, G, Ms, Mt, Ml, Mm, L, m = symbols('x y px py xt yt xm ym xl yl G Ms Mt Ml Mm L m')



h = 1/2 * (px**2 + py**2) - G*Ms/L * 1/(sqrt(x**2 + y**2)) - G*Mt/L * 1/(sqrt((x - xt)**2 + (y - yt)**2)) - G*Ml/L * 1/(sqrt((x - xm)**2 + (y - ym)**2)) - G*Mm/L * 1/(sqrt((x - xl)**2 + (y - yl)**2))

print('Hamiltoniano h:')
print(h)
print('')




dh_dx = -1/L * diff(h, x)
dh_dy = -1/L * diff(h, y)
dh_dpx = 1/m * diff(h, px)
dh_dpy = 1/m * diff(h, py)

print('Derviadas parciales de h:')
print('dh/dx: px_punto =', dh_dx)
print('dh/dy: py_punto =', dh_dy)
print('dh/dpx: vx = ', dh_dpx)
print('dh/dpy: vy = ', dh_dpy)

# Aqui se verá el aporte de energía de cada uno
dh_dxt = diff(h, xt)
dh_dyt = diff(h, yt)

print('Aporte de energía de la Tierra:')
print('dh/dxt:', dh_dxt)
print('dh/dyt:', dh_dyt)
