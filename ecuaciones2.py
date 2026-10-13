from sympy import symbols, diff, sqrt

x, y, px, py, xt, yt, xm, ym, xl, yl, Ms, Mt, Ml, Mm = symbols('x y px py xt yt xm ym xl yl Ms Mt Ml Mm')

funcion = Ms/(sqrt(x**2 + y**2)) + Mt/(sqrt((x - xt)**2 + (y - yt)**2)) + Ml/(sqrt((x - xl)**2 + (y - yl)**2)) + Mm/(sqrt((x - xm)**2 + (y - ym)**2))

dx = px

dy = py

dpx = diff(funcion, x)

dpy = diff(funcion, y)

print('FUNCIONES DE MOVIMIENTO:')
print('dx/dt =', dx)
print('dy/dt =', dy)
print('dpx/dt =', dpx)
print('dpy/dt =', dpy)


