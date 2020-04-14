
![](./orientpy/examples/picture/OrientPy_logo.png)

## Seismic station orientation tools 

Toolbox to determine seismometer orientation using automated (and manual) 
processing of earthquake data. These methods are particularly useful for broadband 
ocean-bottom seismic stations, but are also applicable to broadband land stations
or shorter period instruments (depending on the method selected). Currently the toolbox
includes the following methods: 

- **DL** (Doran and Laske, 2017): Based on Rayleigh-wave polarization at a range of
  frequencies and for the two fundamental mode Rayleigh wave orbits. 

- **BNG** (Braunmiller, Nabelek and Ghods, 2020): Based on P-wave polarization of 
  earthquakes at regional teleseismic distances.  

- **LKSS** (Lim et al., 2018): Based on a harmonic decomposition of radial and 
  transverse receiver functions near zero lag times.

Each method can be used independently to produce an estimate of station orientation, in
terms of the azimuth of seismic component '1' (or 'N').

#### References

- Braunmiller, J., Nabelek, J., and Ghods, A., 2020, Sensor orientation of Iranian broadband
  seismic stations from P-wave particle motion, *Seism. Res. Lett.*, doi:10.1785/0220200019.

- Doran, A. K., and Laske, G., 2017, Ocean-bottom seismometer instrument orientation 
  via automated Rayleigh-wave arrival-angle measurements, *Bull. Seism. Soc. Am.*,
  107, 691-708.

- Lim, H., Kim, Y., Song, T.-R. A., and Shen, Z., 2018, Measurement of 
  seismometer orientation using the tangential P-wave
  receiver function based on harmonic decomposition, *Geophys. J. Int.*, 212,
  1747-1765.

