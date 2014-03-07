from ds9 import *
from scipy import array, exp, ones, optimize, where, median
import pyfits
import operator
from gaussfitter import gaussfit
from pylab import *

# Median Absolute Deviation clipping for input list of numbers.
def mad_clipping(input_data, sigma_clip_level):
    medval = median(input_data)
    sigma = 1.48 * median(abs(medval - input_data))
    high_sigma_clip_limit = medval + sigma_clip_level * sigma
    low_sigma_clip_limit = medval - sigma_clip_level * sigma
    clipped_data = []
    for value in input_data:
        if (value > low_sigma_clip_limit) and (value < high_sigma_clip_limit):
            clipped_data.append(value)
    clipped_data_array = array(clipped_data)
    new_medval = median(clipped_data_array)
    new_sigma = 1.48 * median(abs(medval - clipped_data_array))
    return clipped_data_array, new_medval, new_sigma


def moffat_residuals(params, z, x, y):
    central_height, center_x, center_y, alpha, beta = params
    r2 = (x - center_x)**2 + (y - center_y)**2
    err = z - (central_height * (1 + (r2/(alpha**2)))**(-beta))
    return err
def fit_moffat(center_x, center_y, z, x, y):
    """Returns (central_height, center_x, center_y, alpha, beta
    the moffat parameters of a 2D distribution found by a fit"""
    input_params = [2200.0, center_x, center_y, 2.0, 2.0]
    fit_params, success = optimize.leastsq(moffat_residuals, input_params, 
        args=(z, x, y), maxfev=10000)
    return fit_params

def star_analysis(target_x, target_y):
    snippet_size = 50 # sub_image half_width

    full_image_data = pyfits.getdata(image_filepath)
    
    if (target_x < snippet_size) or (target_x > full_image_data.shape[1] - snippet_size) or (target_y < snippet_size) or (target_y > full_image_data.shape[0] - snippet_size):
        print "You must select a position at least", snippet_size, "pixels from the edge of the image. Aborting."
        return
    
    sub_image_data = full_image_data[target_y-snippet_size:target_y+snippet_size, target_x-snippet_size:target_x+snippet_size]

    # Run the 2D gaussian fit
    gauss_fit_results = gaussfit(sub_image_data)
    # gauss_fit_results has the format:
    # (height, amplitude, x, y, width_x, width_y, rotation angle), image
    gauss_fit_height = gauss_fit_results[0]
    gauss_fit_amplitude = gauss_fit_results[1]
    gauss_fit_x_center = gauss_fit_results[2]
    gauss_fit_y_center = gauss_fit_results[3]
    gauss_fit_x_width = gauss_fit_results[4]
    gauss_fit_y_width = gauss_fit_results[5]
    gauss_fit_x_fwhm = 2.35482*gauss_fit_x_width
    gauss_fit_y_fwhm = 2.35482*gauss_fit_y_width
    gauss_fit_fwhm = (gauss_fit_x_fwhm + gauss_fit_y_fwhm)/2.
    gauss_fit_ecc = (abs(1 - gauss_fit_x_fwhm/gauss_fit_y_fwhm))

    sub_image_data = sub_image_data - gauss_fit_height

    z = []
    x = []
    y = []
    for i in range(len(sub_image_data)):
        i_pixel = i + 0
        for j in range(len(sub_image_data[i])):
            if sub_image_data[i][j] != 0:
                j_pixel = j + 0
                x.append(j_pixel)
                y.append(i_pixel)
                z.append(sub_image_data[i][j])


    r = ((x - gauss_fit_x_center)**2 + (y - gauss_fit_y_center)**2)**0.5

    fit_curve_r = linspace(0, snippet_size, 400)
    gauss_fit_curve = gauss_fit_height + gauss_fit_amplitude * exp(-1* (fit_curve_r)**2 / (2 * ( (gauss_fit_x_width + gauss_fit_y_width)/2.)**2))



    fit_params = fit_moffat(gauss_fit_x_center, gauss_fit_y_center, z, x, y)
    central_height, center_x, center_y, alpha, beta = fit_params
    fwhm_moffat = 2*alpha*sqrt(2**(1/beta) - 1)

    rad = sqrt( (x-center_x)**2 + (y-center_y)**2 )
    h = where( rad < 3 * fwhm_moffat )
    z_cropped = array(z)[h]
    x_cropped = array(x)[h]
    y_cropped = array(y)[h]
    r_cropped = array(r)[h]

    fit_params = fit_moffat(center_x, center_y, z_cropped, x_cropped, y_cropped)

    central_height, center_x, center_y, alpha, beta = fit_params
    fwhm_moffat = 2*alpha*sqrt(2**(1/beta) - 1)

    r2 = fit_curve_r**2
    moffat_fit_curve = (central_height * (1 + (r2/(alpha**2)))**(-beta))


    h_annulus = where((rad > 10 * fwhm_moffat) * (rad < 12 * fwhm_moffat))
    z_annulus = array(z)[h_annulus]
    r_annulus = array(r)[h_annulus]


    clf()
    scatter(r_cropped, z_cropped, color="red", alpha=0.5)
    # scatter(r_annulus, z_annulus, color="orange", alpha=0.1)
    plot(fit_curve_r, gauss_fit_curve, color="blue", label="Gaussian Fit FWHM=" + str(round(gauss_fit_fwhm, 3)))
    plot(fit_curve_r, moffat_fit_curve, color="green", label="Moffat Fit FWHM=" + str(round(fwhm_moffat, 3)))
#     plot(fit_curve_r, median(z_annulus)*ones(fit_curve_r.size), color="green", label="Sky Background=" + str(round(median(z_annulus), 3)))
    xlim(0, fwhm_moffat*3.06)
    ylim(max(gauss_fit_height+gauss_fit_amplitude, central_height)*-0.1, max(gauss_fit_height+gauss_fit_amplitude, central_height)*1.05)
    xlabel("Radius [px]")
    ylabel("Flux Counts [ADU]")
    legend(loc=1)
    show()
    return

def background_analysis(target_x, target_y):
    snippet_size = 50 # sub_image half_width
    full_image_data = pyfits.getdata(image_filepath)
    if (target_x < snippet_size) or (target_x > full_image_data.shape[1] - snippet_size) or (target_y < snippet_size) or (target_y > full_image_data.shape[0] - snippet_size):
        print "You must select a position at least", snippet_size, "pixels from the edge of the image. Aborting."
        return
    sub_image_data = full_image_data[target_y-snippet_size:target_y+snippet_size, target_x-snippet_size:target_x+snippet_size]
    clipped_data_array, new_medval, new_sigma = mad_clipping(sub_image_data.reshape((1, sub_image_data.size))[0], 5)

    print new_medval, new_sigma
    return


print ds9_targets()
d = ds9()

image_filepath = d.get("file")

captured_data = d.get("imexam any coordinate image")


print image_filepath
print captured_data

key_pressed = captured_data.split()[0]
target_x = int(round(float(captured_data.split()[1])))
target_y = int(round(float(captured_data.split()[2])))

print key_pressed, target_x, target_y

if key_pressed == "s":
    star_analysis(target_x, target_y)
    
elif key_pressed == "b":
    background_analysis(target_x, target_y)

    
else:
    print """Key press not recognized. Usage:
        s   select a star for FWHM
        b   select background sky for background statistics
        """




