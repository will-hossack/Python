import matplotlib.pyplot as plt
import math
import tio as t


def decay(t=0):
    cnt = 0
    while cnt < 500:
        cnt += 1
        t += 0.1
        yield t, math.sin(2*math.pi*t)*math.exp(-t/10)

def main():

    ta = []
    ya = []
    for d in decay():
      t,y = d
      ta.append(t)
      ya.append(y)

    plt.plot(ta,ya)
    plt.show()

main()

    
