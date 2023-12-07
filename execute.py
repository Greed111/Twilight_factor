from tools import *

#### QUESTION 1
# Calculate the final evaluation formula with assumption that D is 56, F is 8.
ans1 = mytool1(56, 8, 50)
ans1.calc_V()

# Plot corrected twilight factor for 9 points with pupil distance of 7. If you
# want to plot, uncomment the next line of code.
# paint_z_7()

# The following is worth noting!
# Plot the conjugate diameter and magnification of the corrected twilight coefficient.
# The important point is that adjust the source code when drawing twilight coefficient
# less than 4, for we default to a coefficient outside of 4 to 7 to stop the program.So
# we annotate the next line of code.If you want to use it, you must delete lines 33
# through 38 of tools.py and cancel the annotation of line 32.
# paint_conjugate()


#### QUESTION 2
# Calculate the final CMOS evaluation formula with assumption that D is 56, F is 8.
ans2 = mytool2(56, 8, 50)
ans2.calc_V()

# Plot the signal-to-noise ratio index and quantum efficiency distribution. If you
# want to plot, uncomment the next line of code.
# ans2.plot_Poisson_SNR()

# Plot CMOS twilight vision factor with assumption that D is 56, F is 8. If you
# want to plot, uncomment the next line of code.
# plot_cmos()