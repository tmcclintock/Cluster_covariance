import numpy as np
import scipy.integrate as integ

Dh = 3000. #Mpc/h - Hubble distance

class distances(object):
    def __init__(self, Omega_m, Omega_l):
        self.Om = Omega_m
        self.Ol = Omega_l

    def Ez(self, z):
        return np.sqrt(self.Om*(1+z)**3 + self.Ol)

    def invEz(self, z):
        return 1./self.Ez(z)

    def comoving(self, z):
        z = np.atleast_1d(z)
        return Dh * np.array([integ.quad(self.invEz, 0, zi)[0] for zi in z])

    def angular_dd(self, z):
        return self.comoving(z)/(1+z)

if __name__ == "__main__":
    z = np.linspace(0,3)

    dist = distances(0.3, 0.7)

    import matplotlib.pyplot as plt

    com = dist.comoving(z)
    add = dist.angular_dd(z)
    plt.plot(z, com)
    plt.plot(z, add)
    plt.show()
