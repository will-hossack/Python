import numpy as np
import matplotlib.pyplot as plt


x1 = np.linspace(0.0, 5.0)
x2 = np.linspace(0.0, 2.0)

y1 = np.cos(2 * np.pi * x1) * np.exp(-x1)
y2 = np.cos(2 * np.pi * x2)

fig = plt.figure()

top = fig.add_subplot(2,1,1)
#plt.subplot(2, 1, 1)
top.plot(x1, y1, 'o-')
top.set_title('The top subplot')
top.set_ylabel('Damped oscillation')

bottom = fig.add_subplot(2,1,2)
#plt.subplot(2, 1, 2)

bottom.plot(x2, y2, '.-')
bottom.set_title('The bottom subplot')
bottom.set_xlabel('time (s)')
bottom.set_ylabel('Undamped')

plt.show()
