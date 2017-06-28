def get_spi_spectrum(spi, teff, logg, met):
    """
    Inputs:
        spi -
        teff (K)
        logg
        met
    Output: Spectrum with units
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

    if teff < 2800.:
        teff = 2800
    if logg < (-0.5):
        logg = (-0.5)

    # Giants
    if (teff >= 2500. and teff <= 3500. and logg <= 4.0 and logg >= -0.5):
        func = spi['Cool Giants']
        flux = func.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
    elif (teff >= 4500. and teff <= 5500. and logg <= 4.0 and logg >= -0.5):
        func = spi['Warm Giants']
        flux = func.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
    elif (teff >= 5500. and teff < 6500. and logg <= 4.0 and logg >= -0.5):
        func1 = spi['Warm Giants']
        func2 = spi['Hot Stars']
        flux1 = func1.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
        flux2 = func2.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)

        t_index = (np.abs(gh_teff_overlap - teff)).argmin()
        weight = gh_weights[t_index]
        flux = (flux1*weight + flux2*(1-weight))
    elif (teff >= 3500. and teff < 4500. and logg <= 4.0 and logg >= -0.5):
        func1 = spi['Cool Giants']
        func2 = spi['Warm Giants']
        flux1 = func1.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
        flux2 = func2.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)

        t_index = (np.abs(gc_teff_overlap - teff)).argmin()
        weight = gc_weights[t_index]
        flux = (flux1*weight + flux2*(1-weight))

    # Dwarfs
    elif (teff >= 5500. and teff < 6000. and logg > 4.0):
        func = spi['Warm Dwarfs']
        flux = func.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
    elif (teff >= 2500. and teff <= 3000. and logg > 4.0):
        func = spi['Cool Dwarfs']
        flux = func.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
    elif (teff >= 3000. and teff <= 5500. and logg > 4.0):
        func1 = spi['Cool Dwarfs']
        func2 = spi['Warm Dwarfs']
        flux1 = func1.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
        flux2 = func2.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)

        t_index = (np.abs(d_teff_overlap - teff)).argmin()
        weight = d_weights[t_index]
        flux = (flux1*weight + flux2*(1-weight))

    # Hot stars, have to split this up bcuz of warm stars
    elif (teff >= 6500. and teff <= 12e3 and logg <= 4.0 and logg >= -0.5):
        func = spi['Hot Stars']
        flux = func.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
    elif (teff >= 6000. and teff <= 12e3 and logg > 4.0):
        func = spi['Hot Stars']
        flux = func.get_star_spectrum(logt=np.log10(teff),
                logg=logg, feh=feh)
    else:
        error = ('Parameter out of bounds:'
                 'teff = {0},  logg {1}')
        raise ValueError(error.format(teff, logg))

    return flux
