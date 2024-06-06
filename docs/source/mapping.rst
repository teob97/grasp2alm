Conversion to HEALPix map
=========================

To achieve the conversion into spherical harmonics it is necessary to project the beam into a `HELAPix map <https://en.wikipedia.org/wiki/HEALPix>`_, this allows grasp2alm to use built-in functions such as the `spherical harmonics transform <https://healpy.readthedocs.io/en/latest/healpy_spht.html>`_.

HELAPix map produces a subdivision of a spherical surface in which each pixel covers the same surface area as every other pixel. The resolution is controlled by the parameter :math:`N_{\mathrm{side}}` related to the number of pixel by the simple equation :math:`N_{\mathrm{pixel}}=12N_{\mathrm{side}}^2`.

The discretization performed by Grasp is different from the one used in a HEALPix map, to avoid gaps in the map it is necessary to make an interpolation of the beam at a given theta and phi.