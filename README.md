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

## Attribution 

Please cite [Villaume et al. 2017a](https://arxiv.org/pdf/1705.08906.pdf) if you find this code useful for your research. The BibTex entry for this paper is:

```
@ARTICLE{2017ApJS..230...23V,
       author = {{Villaume}, Alexa and {Conroy}, Charlie and {Johnson}, Benjamin and
         {Rayner}, John and {Mann}, Andrew W. and {van Dokkum}, Pieter},
        title = "{The Extended IRTF Spectral Library: Expanded Coverage in Metallicity, Temperature, and Surface Gravity}",
      journal = {\apjs},
     keywords = {atlases, galaxies: stellar content, infrared: stars, planets and satellites: fundamental parameters, stars: fundamental parameters, techniques: spectroscopic, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Astrophysics of Galaxies},
         year = "2017",
        month = "Jun",
       volume = {230},
       number = {2},
          eid = {23},
        pages = {23},
          doi = {10.3847/1538-4365/aa72ed},
archivePrefix = {arXiv},
       eprint = {1705.08906},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2017ApJS..230...23V},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```


