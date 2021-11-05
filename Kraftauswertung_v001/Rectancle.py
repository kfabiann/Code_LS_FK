import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import RectangleSelector

xxdata = np.linspace(0, 9 * np.pi, num=301)
yydata = np.sin(xxdata) * np.cos(xxdata * 2.4)

a = input("Hallo")

fig, ax = plt.subplots()
line, = ax.plot(xxdata, yydata)
point_max, = ax.plot([], [], marker="o", color="crimson")
point_min, = ax.plot([], [], marker="o", color="crimson")
text_max = ax.text(0, 0, "")
text_min = ax.text(0, 0, "")
xmax = [0.0, 0.0]
xmin = [0.0, 0.0]


def line_select_callback(eclick, erelease):
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata

    mask = (xxdata > min(x1, x2)) & (xxdata < max(x1, x2)) & \
        (yydata > min(y1, y2)) & (yydata < max(y1, y2))
    xmasked = xxdata[mask]
    ymasked = yydata[mask]

    print(xmasked)
    print(ymasked)

    print(len(xmasked))

    if len(xmasked) > 0:
        # xmax = xmasked[np.argmax(ymasked)]
        # ymax = ymasked.max()
        xmax[0] = xmasked.max()
        xmax[1] = ymasked[np.argmax(xmasked)]

        xmin[0] = xmasked.min()
        xmin[1] = ymasked[np.argmin(xmasked)]

        #print(xmax, ymax, xmin, ymin)

        t_max = "xmax: {:.3f}, y(xmax) {:.3f}".format(xmax[0], xmax[1])
        t_min = "xmin: {:.3f}, y(ymin) {:.3f}".format(xmin[0], xmin[1])

        point_min.set_data((xmin[0], xmin[1]))
        point_max.set_data((xmax[0], xmax[1]))

        text_max.set_text(t_max)
        text_max.set_position(((xmax[0] + 0.25), xmax[1]))
        text_min.set_text(t_min)
        text_min.set_position(((xmin[0] - 5), xmin[1]))

        #print(xmax[0], xmax[1])

        fig.canvas.draw_idle()


rs = RectangleSelector(ax, line_select_callback,
                       drawtype='box', useblit=False, button=[1],
                       minspanx=5, minspany=5, spancoords='pixels',
                       interactive=True)

#plt.plot([xmin, ymin], marker="o")


plt.show()

print(xmax[0], xmax[1])
print(xmin[0], xmin[1])
