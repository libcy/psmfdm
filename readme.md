#psmfdm

A cuda based P-SV wave propagation program using hybrid PSM/FDM method. It calculates the horizonal spacial derivatives with pseudo-spectral method, and vertical spacial derivatives with 4th order staggered grid finite difference method. The original program was written in fortran by Yanbin Wang. I rewrote the program in C and then converted it to cuda C, so there's a direct comparison between its CPU and GPU version that illustrates the benefits of GPU. This work was published on [Acta Seismologica Sinica](https://www.researchgate.net/publication/282321018_GPU-based_simulation_of_seismic_wave_propagation_with_hybrid_PSMFDM_method).

* Speedup compared to CPU version(Intel core i5-4570, Nvidia GTX 750)
  ![](https://raw.githubusercontent.com/libcy/psmfdm/master/img/speedup.png)

* Synthetic seismograms in Lanzhou basin model
  ![](https://raw.githubusercontent.com/libcy/psmfdm/master/img/seismogram.png)

* Wavefied snapshots
  ![](https://raw.githubusercontent.com/libcy/psmfdm/master/img/seismogram)
