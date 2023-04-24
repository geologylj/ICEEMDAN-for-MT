# Denoising Application of Magnetotelluric Low-Frequency Signal Processing
This is the origin MATLAB implementation of ICEEMDAN in the following paper:
[Denoising Application of Magnetotelluric Low-Frequency Signal Processing](https://ieeexplore.ieee.org/abstract/document/9904938)

##function
```
	- Run Main.m to see the example
	- acoo.m is aco algorithm
	- addsj.m add triangular wave
	- ceemdan.m is used CEEMDAN to decomposition signal
	- cfdian.m add charge discharge triangular wave
	- Code.m encodes variables into chromosomes for randomly initializing a population
	- colored_noise calculates a colored noise
	- Cross.m completes the crossover operation
	- Decode.m decodes chromosomes
	- EEMD.m is used EEMD to decomposition signal
	- EMD.m computes Empirical Mode Decomposition
	- iceemdan.m is used ICEEMDAN to decomposition signal
```

	
## <span id="citelink">Citation</span>
If you find this repository useful in your research, please consider citing the following paper:

```
@inproceedings{ICEEMDAN-MT-2022,
    author    = {Jin Li and
               Fanhong Ma and
               Jingtian Tang and
               Yecheng Liu and
               Jin Cai},
	journal = {IEEE Transactions on Geoscience and Remote Sensing},
	title     = {Denoising Application of Magnetotelluric Low-Frequency Signal Processing},
    volume    = {60},
    year      = {2022},
}
```
## Contact
If you have any questions, feel free to contact Jin Li through Email (geologylj@163.com) or Github issues. Pull requests are highly welcomed!
