from polygon import Poly
from model import Model
import matplotlib.pylab as plt

p = Poly("a3411-gmrt610mhz-cons.reg")
a3411mod = Model("a3411_redsequence_NumberDensity.fits")
hist, bin_edges = a3411mod.histogram(p, 0.17)

plt.bar(bin_edges[:-1], hist, width=50)
plt.xlim(min(bin_edges), max(bin_edges))
plt.ylim(0, 1.1*max(hist))
plt.savefig("a3411mod_hist.pdf")
