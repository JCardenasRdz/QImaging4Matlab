import os
import matplotlib.pyplot as plt
import pylab as pl
import pydicom
import numpy as np
from roipoly import roipoly
from lorentzian import make_lorentzian_func
from scipy.optimize import curve_fit

#cut and paste directory
dir_for_images="./022316_pH Phantoms/2411_phantom_3_2411_022316__E2_P1_2.16.756.5.5.100.3925202933.28906.1456263202.156/"

os.chdir(dir_for_images)
ConstPixelDims=(int(128), int(128), int(51))
array_dicom=np.zeros(ConstPixelDims,dtype=int)
i=0
for filename in os.listdir(dir_for_images):
	temp = dicom.read_file(filename)
	array_dicom[:,:,i]=temp.pixel_array
	i=i+1


# no method file, hard code ppm offsets:
ppm=np.linspace(-6,6,50)  #these are fucked, go out to 8ppm or so

pl.imshow(array_dicom[:,:,1])
MyROI = roipoly(roicolor='r') # draw new ROI in red color
mask = MyROI.getMask(array_dicom[:,:,1])
Z_signal=[]
for i in range(0,50):
	temp_image=array_dicom[:,:,i]
	Z_signal.append(pl.mean(temp_image[mask]))

#normalize
Z_signal=Z_signal/max(Z_signal)

X0=[.5,.1,.1,.1,2.15,.5,.5,0.5,0,1.8,4.2,5.6]
popt, pcov = curve_fit(make_lorentzian_func(4,False), ppm, Z_signal, X0)
simmed_data=make_lorentzian_func(4,True)(ppm,popt)

print("CEST 4.2: ", popt[2])
print("CEST 5.6: ", popt[3])

plt.plot(ppm,Z_signal,ppm,simmed_data)
plt.show()

wait = input("PRESS ENTER TO END")
