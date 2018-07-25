# SPI-Utils

## Example

```
>>> import matplotlib.pyplot as plt
>>> import spigen

>>> spec = spigen.Spectrum()

>>> teff = 3110.1831261806437
>>> logg = 5.1249530724684771
>>> feh = 0.0

>>> spec = spec.from_coefficients(teff, logg, feh)
>>> plt.plot(spec['wave'], spec['flux'])
>>> plt.show()
```
