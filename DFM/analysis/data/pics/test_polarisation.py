import matplotlib.image as mpimg
import numpy as np
import matplotlib.pyplot as plt
import glob

pic_files = glob.glob('frames/*.jpeg')
pic_files = np.sort(pic_files)
pics = []
intensity_r = []
intensity_g = []
intensity_b = []
for pic_file in pic_files:
    pics.append(mpimg.imread(pic_file))
    intensity_r.append(np.mean(pics[-1][:, :313, 0]) * 255)
    intensity_g.append(np.mean(pics[-1][:, :313, 1]) * 255)
    intensity_b.append(np.mean(pics[-1][:, :313, 2]) * 255)

#plt.imshow(pics[-1][:, :500, :])
plt.plot(intensity_r, 'r-') # / max(intensity_r)
plt.plot(intensity_g, 'g-') # / max(intensity_g)
plt.plot(intensity_b, 'b-') # / max(intensity_b)

plt.savefig('rgb_pol.png')



pic_files = glob.glob('frames_I/*.jpeg')
pic_files = np.sort(pic_files)
pics = []
intensity_r = []
intensity_g = []
intensity_b = []
for pic_file in pic_files:
    pics.append(mpimg.imread(pic_file))
    intensity_r.append(np.mean(pics[-1][:, :, 0]))
    intensity_g.append(np.mean(pics[-1][:, :, 1]))
    intensity_b.append(np.mean(pics[-1][:, :, 2]))

plt.clf()
plt.plot(intensity_r, 'r-') # / max(intensity_r)
plt.plot(intensity_g, 'g-') # / max(intensity_g)
plt.plot(intensity_b, 'b-') # / max(intensity_b)
plt.savefig('rgb_int.png')
