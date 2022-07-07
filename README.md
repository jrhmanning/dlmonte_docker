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

#### Container Parameters

Currently the container will prepare a test simulation in the `/run/` directory which will run a test isotherm of Nitrogen in CuBTC (HKUST-1), outputting a graph and data file to the /interface/ directory for file I/O. This takes a little while to run (ca. 1 hour).

```shell
docker run -v dlm_vol:/interface/ dl_monte
```

## Built With

* DL_MONTE 2.07


## Find Us

* [Gitlab](https://gitlab.com/dl_monte)

## License

This project is licensed under the BSD 3-clause License - see the LICENSE.md file for details.

