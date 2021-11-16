# import numpy as np
#
#
# def pop_dyn_nim(x, y, z, w, an, l, ga, d, rh, my, mz, mw, e):
#     dx = -(an * y + l) * x
#     dy = e * an * x * y + ga * w - (d + my) * y
#     dz = d * y - mz * z
#     dw = rh * z - (ga + mw) * w
#
#     return np.array([dx, dy, dz, dw])
#
#
# def RK4(f, I, y0, n):
#     h = (I[1] - I[0])/n
#     t = I[0]
#     y = y0
#     for i in range(n):
#         k1 = h * f(t, y0)
#         k2 = h * f(t + h/2, y0 + k1/2)
#         k3 = h * f(t h/2, y0+ k2/2)
#         k4 = h * f(t + h, y0+k3)
#         y1 = y0 + (k1 + 2 * (k2 + k3) + k4) / 6
#         y = np.array(y0, y1)