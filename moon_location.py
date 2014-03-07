import ephem
from math import pi
from pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

observatory = ephem.Observer()

# Kitt Peak
# observatory.long = "-111.59666667"
# observatory.lat = "31.9583333333"
# observatory.elevation = 2095

# Lick
# observatory.long = "-121.642847"
# observatory.lat = "37.341717"
# observatory.elevation = 1283

# Cerro Tololo
# observatory.long = "-70.804"
# observatory.lat = "-30.169"
# observatory.elevation = 2123

# PAIRITEL
observatory.long = "-110.87904"
observatory.lat = "31.68084"
observatory.elevation = 2616


dblist = []
dbfile = file("rrlyrae_coordinates_2.txt", "r")
for line in dbfile:
    object_name = line.split("\t")[0]
    ra_deg = float(line.split("\t")[1])
    ra_hr = float(line.split("\t")[1]) / 15.0
    dec_deg = float(line.split("\t")[2])
    new_point = ephem.Equatorial(str(ra_hr), str(dec_deg), epoch="2000")
    ephem_line = (object_name + ",f|M|F7," + str(new_point.ra) + "," + 
        str(new_point.dec) + ",0,2000")
    dblist.append(ephem_line)
dbfile.close

sun = ephem.Sun()
moon = ephem.Moon()

final_output_file = file("observing_dates.txt", "w")
final_output_file.write("object\tstart_date\tend_date\twindow_size\tmin_moon_sep\n")
final_output_file.close()

big_list_of_observing_window_sizes = []
big_list_of_adjusted_daily_dates = []

start_date = ephem.date("2011/02/01")
end_date = ephem.date("2011/04/30")

alt_constraint = 40.0 # minimum altitude in degrees defining the "good" window
window_constraint = 6.0 # minimum length in hours of "good" window

for obj in dblist:
    max_alt = 0
    max_alt_date = ephem.date("2000")
    star = ephem.readdb(obj)
    daily_dates = []
    observing_window_sizes = []
    min_moon_spearations = []
    for n in range(int(end_date - start_date)):
        observatory.date = start_date + n*1
        observatory.date = observatory.next_antitransit(sun)
        star.compute(observatory)
        alt = star.alt*180/pi
        if alt > max_alt:
            max_alt = alt
            max_alt_date = observatory.date
        dates = []
        altitudes = []
        good_dates = []
        moon_separations = []
        mid_date = observatory.date
        for m in range(500):
            observatory.date = ephem.date(float(mid_date) - 0.5 + m/500.)
            star.compute(observatory)
            new_alt = float(star.alt*180/pi)
            sun.compute(observatory)
            sun_alt = float(sun.alt*180/pi)
            new_date = float(observatory.date) - float(mid_date)
            dates.append(new_date)
            altitudes.append(new_alt)
            moon.compute(observatory)
            moon_separation = ephem.separation(moon, star)
            if (new_alt > alt_constraint) and (sun_alt < -18):
                good_dates.append(new_date)
                moon_separations.append(moon_separation)
        try:
            total_good_time = good_dates[-1] - good_dates[0]
        except(IndexError):
            total_good_time = 0
        daily_dates.append(float(mid_date))
        observing_window_sizes.append(24.*total_good_time)
        try:
            min_moon_spearations.append(57.2957796*float(min(moon_separations)))
        except:
            min_moon_spearations.append(181)
#         plot_title = (star.name + " around " + str(ephem.date(float(mid_date))))
#         fig = plt.figure(figsize=(8, 8))
#         ax1 = fig.add_subplot(1,1,1)
#         ax1.plot(dates, altitudes, 
#             marker = "s", color = "red", linestyle="none", label = "Path")
#         ax1.set_ylim(-90, 90)
#         ax1.set_xlabel("PyEphem Date")
#         ax1.set_ylabel("Altitude")
#         ax1.set_title(plot_title)
#         canvas = FigureCanvas(fig)
#         canvas.print_figure(star.name + "_around_" + str(float(mid_date)) + ".png", dpi=240)
    
    final_daily_dates = []
    for d in range(len(daily_dates)):
        if observing_window_sizes[d] > window_constraint:
            final_daily_dates.append(daily_dates[d])
    adjusted_daily_dates = []
    for date in daily_dates:
        adjusted_daily_dates.append(date - float(start_date))
    
    plot_title = (star.name + " observing window sizes (hours above " + 
        str(alt_constraint) + " deg in altitude)")
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(adjusted_daily_dates, observing_window_sizes, 
        marker = "s", color = "red", linestyle="-", label = "Path")
    ax1.set_ylim(0, 8)
    ax1.set_xlim(0, int(end_date - start_date))
    ax1.set_xlabel("Semester Date")
    ax1.set_ylabel("Window Size [hours]")
    ax1.set_title(plot_title)
    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(adjusted_daily_dates, min_moon_spearations, 
        marker = "o", color = "blue", linestyle="-", label = "Path")
    ax2.set_ylim(0, 180)
    ax2.set_xlim(0, int(end_date - start_date))
    ax2.set_xlabel("Semester Date")
    ax2.set_ylabel("Minimum Moon Separation [deg]")
    canvas = FigureCanvas(fig)
    canvas.print_figure("plots/" + star.name + "_observing_windows.png", dpi=240)
    
#     plot_title = (star.name + " minimum nightly moon spearations")
#     fig = plt.figure(figsize=(8, 8))
#     ax1 = fig.add_subplot(1,1,1)
#     ax1.plot(adjusted_daily_dates, min_moon_spearations, 
#         marker = "s", color = "blue", linestyle="-", label = "Path")
#     ax1.set_ylim(0, 200)
#     ax1.set_xlim(0, int(end_date - start_date))
#     ax1.set_xlabel("Semester Date")
#     ax1.set_ylabel("Minimum Moon Separation [deg]")
#     ax1.set_title(plot_title)
#     canvas = FigureCanvas(fig)
#     canvas.print_figure("plots/" + star.name + "_moon_separations.png", dpi=240)

    big_list_of_observing_window_sizes.append(observing_window_sizes)
    big_list_of_adjusted_daily_dates.append(adjusted_daily_dates)

    print star.name, max_alt_date, max_alt
    try:
        print ("Observing windows of > " + str(window_constraint) + 
            " hours at > " + str(alt_constraint) + " deg altitude occur " + 
            "from %s to %s for a total of %f days." % (str(ephem.date(final_daily_dates[0])),
            str(ephem.date(final_daily_dates[-1])), final_daily_dates[-1] - final_daily_dates[0]))
        print ("Minimum moon separation during tested dates is " + 
            str(min(min_moon_spearations)) + " deg.")
        final_output_file = file("observing_dates.txt", "a")
        final_output_file.write(str(star.name) + "\t" + 
            str(ephem.date(final_daily_dates[0])).split()[0] + "\t" + 
            str(ephem.date(final_daily_dates[-1])).split()[0] + "\t" + 
            str(int(final_daily_dates[-1] - final_daily_dates[0])) + "\t" +
            str(min(min_moon_spearations)) + "\n")
        final_output_file.close()
    except:
        print("No good windows :(")
        final_output_file = file("observing_dates.txt", "a")
        final_output_file.write(str(star.name) + "\t" + 
            str(0) + "\t" + 
            str(0) + "\t" + 
            str(0) + "\n")
        final_output_file.close()

plot_title = ("Observing window sizes (hours above " + str(alt_constraint) 
    + " deg in altitude)")
fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(1,1,1)
for n in range(len(big_list_of_adjusted_daily_dates)):
    ax1.plot(big_list_of_adjusted_daily_dates[n], 
        big_list_of_observing_window_sizes[n], 
        marker = "", color = "red", linestyle="-", label = "Path")
ax1.set_ylim(0, 8)
ax1.set_xlim(0, int(end_date - start_date))
ax1.set_xlabel("Semester Date")
ax1.set_ylabel("Window Size [hours]")
ax1.set_title(plot_title)
canvas = FigureCanvas(fig)
canvas.print_figure("observing_windows.png", dpi=240)

for n in range(len(big_list_of_adjusted_daily_dates[0])):
    day = n
    num_good_targets = 0
    for m in range(len(big_list_of_observing_window_sizes)):
        if big_list_of_observing_window_sizes[m][n] > window_constraint:
            num_good_targets += 1
    print "On day " + str(ephem.date(day + start_date)) + " the number of good targets is " + str(num_good_targets)