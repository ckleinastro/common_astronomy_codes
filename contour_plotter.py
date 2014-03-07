from numpy import *
from scipy.interpolate import interp2d
import matplotlib.cm as cm
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import mpl_toolkits.mplot3d.axes3d as p3
import sys
import scipy.stats as stats

def percent_level(z, low_val):
    n_pix_total = float(z.shape[0]*z.shape[1])
    n_pix_above = float(where(z>low_val)[0].shape[0])
    return float(n_pix_above/n_pix_total)

def find_level(z, target_percent):
    n_pix_total = float(z.shape[0]*z.shape[1])
    new_level = z.mean()/(n_pix_total)
    current_percent_level = percent_level(z, new_level)
    while abs(current_percent_level - target_percent) > 0.001:
        if current_percent_level < target_percent:
            new_level = new_level*0.99
        else:
            new_level = new_level*1.01
        current_percent_level = percent_level(z, new_level)
    return new_level

data = loadtxt("ogle_v_per_amp.txt", delimiter="|")
x_data = array([data[:,0]]).T
y_data = array([data[:,1]]).T

x_min = 0.2
x_max = 1.0
y_min = 0.0
y_max = 1.6

# Create some dummy data
rvs = np.append(x_data, y_data, axis=1)

kde = stats.kde.gaussian_kde(rvs.T)

# Regular grid to evaluate kde upon
x_flat = np.r_[x_min:x_max:256j]
y_flat = np.r_[y_min:y_max:256j]
x,y = np.meshgrid(x_flat,y_flat)
grid_coords = np.append(x.reshape(-1,1),y.reshape(-1,1),axis=1)

z = kde(grid_coords.T)
z = z.reshape(256,256)



l_999 = find_level(z, 0.001)
l_99 = find_level(z, 0.01)
l_975 = find_level(z, 0.025)
l_95 = find_level(z, 0.05)
l_90 = find_level(z, 0.10)
l_85 = find_level(z, 0.15)
l_80 = find_level(z, 0.20)
l_70 = find_level(z, 0.30)
levels_list=[l_70, l_80, l_85, l_90, l_95, l_975, l_99, l_999, z.max()]

fig = plt.figure(figsize=(8.265, 8.0))
ax1 = fig.add_subplot(1,1,1)
plot_title = ("Contour Plot")
# ax1.contourf(x, y, z, levels=levels_list, origin="lower", cmap=cm.gist_yarg)
ax1.contour(x, y, z, levels=levels_list, origin="lower", colors='k')
ax1.scatter(x_data, y_data, alpha=0.15, color="blue", linewidths=0)
ax1.set_xlabel("Period")
ax1.set_ylabel("Amplitude")
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_title(plot_title)
canvas = FigureCanvas(fig)
canvas.print_figure("contour.png", dpi=300)
close("all")
