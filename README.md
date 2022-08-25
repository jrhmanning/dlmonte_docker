# dlmonte_docker
Mirror of the docker hub image https://hub.docker.com/repository/docker/jrhmanning/dl_monte

# DL_MONTE Docker

This container houses the executable code for [DL_MONTE 2.07](https://gitlab.com/dl_monte), a general-purpose Monte Carlo simulation program for a range of molecular simulations.  Further information on DL_MONTE can be found [here](https://gitlab.com/dl_monte/user-hub/-/wikis/home).

## Getting Started

These instructions will cover usage information and for the docker container 

### Prerequisites

In order to run this container you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

The container has several pre-formatted template simulation files and runscripts contained within `/scripts/`, which can be run directly when calling a script with the container:

```shell
docker run -v $PWD/interface/:/run/interface/ jrhmanning/dl_monte:latest python isotherm_runner.py --InputFolder=/run/interface --OutputFolder=/run/interface --Framework=Cu_BTC
```

#### Simulation templates
Currently, the repository contains a single template control object:

- `AdsorptionExample`, within the `isotherm_control_generator.py` file. 
  - This describes a simple GCMC isotherm with no advanced sampling, which can be modified to suit most standard adsorption simulation workflows. 

#### Simulation scripts
The repository contains the following simulation run scripts, which uses `dlmontepython.simtask` to automate GCMC tasks:

- `isotherm_runner.py`
  - Runs an isotherm on a fixed-atom framework found in a `.cif` file. Runs through different pressure values to perform GCMC, and return a `.csv` and `.png` file summarising the results.

#### Simulation Parameters

Simulation parameters can be controlled using argument parsing inside the various python scripts to control or amend parameters. These should be uniform across all simulation scripts included.

* `InputFolder`
  * The input directory containing your input `.cif` file
* `OutputFolder`
  * The output directory whee your simulation outputs will be located. 
* `FrameworkName`
  * The name of your `.cif` file (without the file type, e.g. `Cu_BTC`, not `Cu_BTC.cif`)
* `GasComposition`
  * The name and relative frequency of each proble molecule (N.B. only single gas molecules is currently supported)
* `Temperature`
  * The single temperature of your isotherm (in K)
* `Pressures`
  * The pressure values if your isotherm, as a comma-separated string (e.g. `'1,2,5,1000'`)
* `Charges`
  * A boolean to turn off the Ewald summation, if you want to run a much faster simulation


## Built With

* DL_MONTE 2.07
* dlmontepython
* [ase](https://wiki.fysik.dtu.dk/ase/)


## Find Us

* [Gitlab](https://gitlab.com/dl_monte)

## License

This project is licensed under the BSD 3-clause License - see the LICENSE.md file for details.

