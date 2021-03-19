# stress-domains

Construct a deformation mechanism map for quartz with the following constitutive equations (from Passchier and Trouw, 1996).  The x-axis should be homologous temperatures from 0.2 to 1.0 and the y-axis log10($\sigma$/$\mu$) fromm 0 to -6.  Plot strain rate contours 10^-6 through 10^-15. 

## Usage

### Set Up Environment

### Clone Repository

The following will create a directory called `stress-domains` in the current working directory.

```shell
$ git clone https://github.com/acbuie/stress-domains.git
```

### Create Virtual Environment

Navigate to directory where you cloned the repository. 

Mac: 
```shell
$ python3 -m venv .venv
$ source .venv/bin/activate
$ pip3 install -r requirements.txt
```

For more information, see: https://docs.python.org/3/library/venv.html

## Creating the Plot

Edit the plot config with: 

```shell
$ nano parameters.py
```

Default values are set in `parameters.py`. 

To run:

```shell
$ python3 main.py
```