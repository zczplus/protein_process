import numpy as np
import cv2 as cv
import matplotlib.pyplot as plt

if __name__ == '__main__':
    img = cv.imread('img/img1_TE9_255/TE9.jpg', 0)
    hist, bins = np.histogram(img.flatten(), 256, [0, 256])

    cdf = hist.cumsum()

    cdf_m = np.ma.masked_equal(cdf, 0)
    cdf_m = (cdf_m - cdf_m.min()) * 255 / (cdf_m.max() - cdf_m.min())
    cdf = np.ma.filled(cdf_m, 0).astype('uint8')

    img2 = cdf[img]
    hist, bins = np.histogram(img2.flatten(), 256, [0, 256])
    cdf = hist.cumsum()
    cdf_normalized = cdf * float(hist.max()) / cdf.max()

    plt.plot(cdf_normalized, color='b')
    plt.hist(img2.flatten(), 256, [0, 256], color='r')
    plt.xlim([0, 256])
    plt.legend(('cdf', 'histogram'), loc='upper left')
    plt.show()

    pass
