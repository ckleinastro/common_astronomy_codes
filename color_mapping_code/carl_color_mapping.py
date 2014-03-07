from scipy import loadtxt, zeros, arange
from scipy.interpolate import interp1d
from pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.ticker import MultipleLocator


# This data was extracted from Carl's plot of his colormap
raw_red_color = loadtxt("red_mapping.txt")
raw_green_color = loadtxt("green_mapping.txt")
raw_blue_color = loadtxt("blue_mapping.txt")

red_function = interp1d(raw_red_color[:,0], raw_red_color[:,1], kind="linear")
green_function = interp1d(raw_green_color[:,0], raw_green_color[:,1], kind="linear")
blue_function = interp1d(raw_blue_color[:,0], raw_blue_color[:,1], kind="linear")


def make_color(val, brightness=1.0):
    # "val" is the value that we want to map to color space
    # "brightness" is how bright we want the color to be
    val *= 255
    color = [red_function(val)/255. * brightness, green_function(val)/255. * brightness, blue_function(val)/255. * brightness, 1.0]
    return color

data = loadtxt("W1_plot_data.txt")
# W1 data for RRL stars (143)
# period, absolute_mag, metallicity, amplitude

data_to_colorize = data[:,2]    # colorize by metallicity
data_to_darken = data[:,3]      # darken/brighten by amplitude

color_vals = (data_to_colorize - data_to_colorize.min()) / (data_to_colorize.max() - data_to_colorize.min())
brightness_vals = (data_to_darken - data_to_darken.min()) / (data_to_darken.max() - data_to_darken.min())

color_list=[]
for n in range(len(color_vals)):
    color_list.append(make_color(color_vals[n], brightness_vals[n]))

# blackness is a 256x256 array with 1-0 going down the columns
blackness = zeros((256, 256))
for i in range(256):
    blackness[:,i] = linspace(0, 255, 256)/255.
# now create the 2D colormap, multiplying by blackness to dimm vertically
color_display = zeros((256, 256, 3))
color_display[:,:,0] = (red_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.) * blackness
color_display[:,:,1] = (green_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.) * blackness
color_display[:,:,2] = (blue_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.) * blackness


close("all")
fig = plt.figure(figsize=(8, 8), dpi=144)
fig.suptitle("W1 Period-Luminosity Relation with Metallicity and Amplitude Color Mapping")
    
ax1 = fig.add_subplot(2,1,2)
ax1.scatter(log10(data[:,0]/data[:,0].mean()), data[:,1], c=color_list, marker='o', s=50, linewidth=0, alpha=1.0)
data_range = data[:,1].max() - data[:,1].min()
ax1.set_ylim(data[:,1].max() + 0.1*data_range, data[:,1].min() - 0.1*data_range)

data_width = log10(data[:,0]/data[:,0].mean()).max() - log10(data[:,0]/data[:,0].mean()).min()
ax1.set_xlim(log10(data[:,0]/data[:,0].mean()).min() - 0.1*data_width, log10(data[:,0]/data[:,0].mean()).max() + 0.1*data_width)

ax1.set_ylabel("W1 Absolute Magnitude")
ax1.set_xlabel("Log(Period/P_0)")

ax2 = fig.add_subplot(2,1,1)
im = ax2.imshow(color_display, origin="lower", interpolation="lanczos", 
    extent=[data_to_colorize.min(), data_to_colorize.max(), data_to_darken.min(), data_to_darken.max()])
ax2.set_xlabel("Metallicity")
ax2.xaxis.set_label_position("top")
ax2.set_ylabel("Amplitude")

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.5)
minorLocator_x = MultipleLocator(0.05)
ax2.xaxis.set_major_locator(majorLocator_x)
ax2.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.1)
minorLocator_y1 = MultipleLocator(0.05)
ax2.yaxis.set_major_locator(majorLocator_y1)
ax2.yaxis.set_minor_locator(minorLocator_y1)

# This allows us to set the font size of the tick number labels. 
x_gridlines = ax2.xaxis.get_gridlines()
for tl in ax2.get_xticklabels():
    tl.set_fontsize(12)

# set the tick labels on for both the top&bottom (x) and left&right (y)
for tick in ax2.yaxis.get_major_ticks():
    tick.label1On = True
    tick.label2On = True
for tick in ax2.xaxis.get_major_ticks():
    tick.label1On = False
    tick.label2On = True

# pos = [left, bottom, width, height]
ax2.set_position([0.1, 0.79, 0.8, 0.1])
ax1.set_position([0.1, 0.07, 0.8, 0.7])


show()
canvas = FigureCanvas(fig)
canvas.print_figure("colormapping.pdf", dpi=144)
close("all")