import os
import sys
from pkg_resources import resource_filename
import numpy as np

class Spectrum(object):

    def __init__(self):
        """
        Instantiate the Spectrum class by reading in the
        coefficients
        """
        cd_file = resource_filename('spigen', 'Coefficients/Cool_Dwarfs.dat')
        cg_file = resource_filename('spigen', 'Coefficients/Cool_Giants.dat')
        wd_file = resource_filename('spigen', 'Coefficients/Warm_Dwarfs.dat')
        wg_file = resource_filename('spigen', 'Coefficients/Warm_Giants.dat')
        hs_file = resource_filename('spigen', 'Coefficients/Hot_Stars.dat')

        cd = np.genfromtxt(cd_file)
        cg = np.genfromtxt(cg_file)
        wd = np.genfromtxt(wd_file)
        wg = np.genfromtxt(wg_file)
        hs = np.genfromtxt(hs_file)

        cd, cg, wd, wg, hs = cd.T, cg.T, wd.T, wg.T, hs.T

        # Wavelength is the first point, remove it
        # but first save it as an attribute
        self.wave = cd[0,:]

        self.cd = cd[1:len(cd),:]
        self.cg = cg[1:len(cg),:]
        self.wd = wd[1:len(wd),:]
        self.wg = wg[1:len(wg),:]
        self.hs = hs[1:len(hs),:]

    def cool_dwarfs(self, feh, logt, logg):
        log_flux = (self.cd[1,:] +
                    self.cd[2,:]*logt +
                    self.cd[3,:]*feh +
                    self.cd[4,:]*logg +
                    self.cd[5,:]*feh*feh +
                    self.cd[6,:]*logt*logt +
                    self.cd[7,:]*logg*logg +
                    self.cd[8,:]*logt*feh +
                    self.cd[9,:]*logt*logg +
                    self.cd[10,:]*feh*logg +
                    self.cd[11,:]*feh*feh*feh +
                    self.cd[12,:]*logt*logt*logt +
                    self.cd[13,:]*logg*logg*logg +
                    self.cd[14,:]*logt*logt*feh +
                    self.cd[15,:]*logt*feh*feh +
                    self.cd[16,:]*logg*logt*logt +
                    self.cd[17,:]*logt*logt*logt*logt +
                    self.cd[18,:]*feh*feh*feh*feh +
                    self.cd[19,:]*logt*logt*feh*feh +
                    self.cd[20,:]*logt*logt*logt*feh +
                    self.cd[21,:]*logt*logt*logt*logt*logt
                    )

        return self.cd[0,:] * np.exp(log_flux)

    def cool_giants(self, feh, logt, logg):
        log_flux = (self.cg[1,:] +
                    self.cg[2,:]*logt +
                    self.cg[3,:]*feh +
                    self.cg[4,:]*logg +
                    self.cg[5,:]*logt*logt +
                    self.cg[6,:]*logg*logg +
                    self.cg[7,:]*feh*feh +
                    self.cg[8,:]*feh*logg +
                    self.cg[9,:]*logt*logg +
                    self.cg[10,:]*logt*feh +
                    self.cg[11,:]*logt*logt*logt +
                    self.cg[12,:]*logg*logg*logg +
                    self.cg[13,:]*feh*feh*feh +
                    self.cg[14,:]*logt*logg*feh +
                    self.cg[15,:]*logt*logt*feh +
                    self.cg[16,:]*logt*logt*logg +
                    self.cg[17,:]*feh*feh*logt +
                    self.cg[18,:]*feh*feh*logg +
                    self.cg[19,:]*logt*logg*logg +
                    self.cg[20,:]*feh*logg*logg +
                    self.cg[21,:]*logt*logt*logt*logt
                    )

        return self.cg[0,:] * np.exp(log_flux)

    def warm_giants(self, feh, logt, logg):
        log_flux = (self.wg[1,:] +
                    self.wg[2,:]*logt +
                    self.wg[3,:]*feh +
                    self.wg[4,:]*logg +
                    self.wg[5,:]*logt*logt +
                    self.wg[6,:]*logg*logg +
                    self.wg[7,:]*feh*feh +
                    self.wg[8,:]*logt*feh +
                    self.wg[9,:]*logt*logg +
                    self.wg[10,:]*logg*feh +
                    self.wg[11,:]*logt*logt*logt +
                    self.wg[12,:]*logg*logg*logg +
                    self.wg[13,:]*feh*feh*feh +
                    self.wg[14,:]*logt*logt*feh +
                    self.wg[15,:]*logt*feh*feh +
                    self.wg[16,:]*logg*logt*logt +
                    self.wg[17,:]*logg*logg*logt +
                    self.wg[18,:]*logt*logt*logt*logt +
                    self.wg[19,:]*feh*feh*feh*feh +
                    self.wg[20,:]*logt*logt*feh*feh +
                    self.wg[21,:]*logt*logt*logg*logg +
                    self.wg[22,:]*feh*feh*logg*logg +
                    self.wg[23,:]*logt*logt*logt*logt*logt
                    )

        return self.wg[0,:] * np.exp(log_flux)

    def warm_dwarfs(self, feh, logt, logg):
        log_flux = (self.wd[1,:] +
                    self.wd[2,:]*logt +
                    self.wd[3,:]*feh +
                    self.wd[4,:]*logg +
                    self.wd[5,:]*logt*logt +
                    self.wd[6,:]*logg*logg +
                    self.wd[7,:]*feh*feh +
                    self.wd[8,:]*logt*feh +
                    self.wd[9,:]*logt*logg +
                    self.wd[10,:]*logt*logt*logt +
                    self.wd[11,:]*logt*logg*logg +
                    self.wd[12,:]*feh*feh*feh +
                    self.wd[13,:]*logt*logt*logg +
                    self.wd[14,:]*logt*logt*feh +
                    self.wd[15,:]*logt*feh*feh +
                    self.wd[16,:]*logt*logg*feh +
                    self.wd[17,:]*logg*feh*feh +
                    self.wd[18,:]*logt*logt*logt*logt +
                    self.wd[19,:]*logg*logg*logg*logg +
                    self.wd[20,:]*logt*logt*logt*logg +
                    self.wd[21,:]*feh*logt*logt*logt +
                    self.wd[22,:]*feh*feh*logt*logt +
                    self.wd[23,:]*feh*feh*feh*logt +
                    self.wd[24,:]*logt*logt*logg*logg +
                    self.wd[25,:]*feh*logt*logt*logg +
                    self.wd[26,:]*logt*logt*logt*logt*logt
                    )

        return self.wd[0,:] * np.exp(log_flux)

    def hot_stars(self, feh, logt, logg):
        log_flux = (self.hs[1,:] +
                    self.hs[2,:]*logt +
                    self.hs[3,:]*feh +
                    self.hs[4,:]*logg +
                    self.hs[5,:]*logt*logt +
                    self.hs[6,:]*feh*feh +
                    self.hs[7,:]*logg*logg +
                    self.hs[8,:]*logt*logg +
                    self.hs[9,:]*logt*feh +
                    self.hs[10,:]*logg*feh +
                    self.hs[11,:]*logt*logt*logt +
                    self.hs[12,:]*logg*logg*logg +
                    self.hs[13,:]*feh*feh*feh +
                    self.hs[14,:]*logt*logg*feh +
                    self.hs[15,:]*logt*logt*feh +
                    self.hs[16,:]*logt*logt*logg +
                    self.hs[17,:]*feh*feh*logt +
                    self.hs[18,:]*feh*feh*logg +
                    self.hs[19,:]*logt*logg*logg +
                    self.hs[20,:]*feh*logg*logg +
                    self.hs[21,:]*logt*logt*logt*logt
                    )

        return self.hs[0,:] * np.exp(log_flux)

    def from_coefficients(self, teff, logg, feh):
        """ Generate a stellar spectrum using the
        coeffecients

        Parameters:
        -----------
        teff: float
            effective temperture
        logg: float
            surface gravity
        feh: float
            metallicity

        Output:
        -------
        flux: ndarray
            interpolated spectrum
        """

        """
        These weights are used later to ensure
        smooth behavior.
        """
        # Overlap of cool dwarf and warm dwarf training sets
        d_teff_overlap = np.linspace(3000, 5500, num=100)
        d_weights = np.linspace(1, 0, num=100)

        # Overlap of warm giant and hot star training sets
        gh_teff_overlap = np.linspace(5500, 6500, num=100)
        gh_weights = np.linspace(1, 0, num=100)

        # Overlap of warm giant and cool giant training sets
        gc_teff_overlap = np.linspace(3500, 4500, num=100)
        gc_weights = np.linspace(1, 0, num=100)

        """
        Setting up some boundaries
        """
        teff2 = teff
        logg2 = logg
        if teff2 <= 2800.:
            teff2 = 2800
        if logg2 < (-0.5):
            logg2 = (-0.5)

        # Normalizing to solar values
        logt = np.log10(teff2) - 3.7617
        logg = logg - 4.44

        # Giants
        if (teff2 >= 2500. and teff2 <= 3500. and logg2 <= 4.0 and logg2 >= -0.5):
            flux = self.cool_giants(feh, logt, logg)

        elif (teff2 >= 4500. and teff2 <= 5500. and logg2 <= 4.0 and logg2 >= -0.5):
            flux = self.warm_giants(feh, logt, logg)

        elif (teff2 >= 5500. and teff2 < 6500. and logg2 <= 4.0 and logg2 >= -0.5):
            flux1 = self.warm_giants(feh, logt, logg)
            flux2 = self.hot_stars(feh, logt, logg)

            t_index = (np.abs(gh_teff_overlap - teff2)).argmin()
            weight = gh_weights[t_index]
            flux = (flux1*weight + flux2*(1-weight))

        elif (teff2 >= 3500. and teff2 < 4500. and logg2 <= 4.0 and logg2 >= -0.5):
            flux1 = self.cool_giants(feh, logt, logg)
            flux2 = self.warm_giants(feh, logt, logg)

            t_index = (np.abs(gc_teff_overlap - teff2)).argmin()
            weight = gc_weights[t_index]
            flux = (flux1*weight + flux2*(1-weight))

        # Dwarfs
        elif (teff2 >= 5500. and teff2 < 6000. and logg2 > 4.0):
            flux = self.warm_dwarfs(feh, logt, logg)

        elif (teff2 >= 2500. and teff2 <= 3000. and logg2 > 4.0):
            flux = self.cool_dwarfs(feh, logt, logg)

        elif (teff2 >= 3000. and teff2 <= 5500. and logg2 > 4.0):
            flux1 = self.cool_dwarfs(feh, logt, logg)
            flux2 = self.warm_dwarfs(feh, logt, logg)

            t_index = (np.abs(d_teff_overlap - teff2)).argmin()
            weight = d_weights[t_index]
            flux = (flux1*weight + flux2*(1-weight))

        # Hot stars, have to split this up bcuz of warm stars
        elif (teff2 >= 6500. and teff2 <= 12e3 and logg2 <= 4.0 and logg2 >= -0.5):
            flux = self.hot_stars(feh, logt, logg)

        elif (teff2 >= 6000. and teff2 <= 12e3 and logg2 > 4.0):
            flux = self.hot_stars(feh, logt, logg)
        else:
            error = ('Parameter out of bounds:'
                     'teff = {0},  logg {1}')
            raise ValueError(error.format(teff2, logg))

        spec = {}
        i = (self.wave >= 0.36)
        spec['wave'] = self.wave[i]
        spec['flux'] = flux[i]

        return spec

if __name__=='__main__':
    pass
