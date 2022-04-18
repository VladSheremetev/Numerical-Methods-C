import numpy as np
import matplotlib.pyplot as plt
x = np.arange(-1.0, 1.0, 0.01)
y = np.arange(-1.0, 1.0, 0.01)
y1 = y**3
y2 = np.tan(x)
plt.plot(y1, y, 'r', label='f1')
plt.plot(x, y2, 'b', label='f2')
plt.plot(0.734, 0.902, 'o')
plt.plot(0.0, 0.0, 'o')
plt.plot(-0.734, -0.902, 'o')
plt.annotate('A(0.734, 0.902)', (0.734 + 0.02, 0.902 - 0.1))
plt.annotate('B(0.0, 0.0)', (0.0 + 0.03, 0.0 - 0.05))
plt.annotate('C(-0.734, -0.902)', (-0.734 + 0.03, -0.902 - 0.07))
plt.xlabel('x')
plt.ylabel('y')
plt.title('y(x)')
plt.grid()
plt.legend()
plt.show()

x = np.arange(0.0, 1.1, 0.01)
y = np.arange(-1.0, 1.0, 0.01)

xg, yg = np.meshgrid(x, y)

zg = (xg - yg**3)**2 + (np.tan(xg) - yg)**2

plt.contour(xg, yg, zg, levels = 100)
plt.title("Level lines of squared discrepancy sum")
plt.colorbar()
plt.show()
