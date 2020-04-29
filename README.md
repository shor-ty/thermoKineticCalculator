The open source thermo-kinetic calculator (TKC) initiated by Tobias Holzmann.

### Author ###
* Tobias Holzmann (inventor and main developer)

### Contributors ###
* None

### The open source thermo-kinetic calculator library###
This open source projects is directed to all people who want to deal with detailed chemical reactions. The repository gives one the freedom to use detailed chemistry in any appliation. The application consists of different libraries which handles *thermodynamics* and *chemistry* data while using *gas kinetic* theory.

Additionally, the repository shares ready to use applications for:
* Transient ideal-homogeneous reactors calculations
* Steady-State flamelet calculator
* Others will follow based on the interested

### Why a new library? ###
If one wants to use thermo-kinetic calculations, there are already a wide range of applications and libraries available such as the *Cantera* and *FlameMaster* projects. However, in the year 2014, Tobias was investigating a lot into the flamelet model using OpenFOAM. For that purpose, only the binary flamelet generator built by the CRECK-Modeling group from Milano were available and Tobias was too lazy to investigate into *Cantera* or *FlameMaster*. However, the main driving force was that Tobias wanted to get more familiar with implementing theory into c++ and also wanted to increase his c++ knowledge during his Ph.D.


### History ###
Since 2014, Tobias was working a lot on the code only during his spare time. Since 2017, the project was not further investigated as time was limited and other things were more important. In 2020, Tobias reinvestigated into the code, changed and re-organized a lot of stuff as well as made a huge simplification.

### Status Quo ###
The thermo-kinetic library consists of different classes that are of up to date and old (depreciated). Right now the only classes are of interest:
* src/definitions
* src/manipulation
* src/thermoKinetics
* src/mathematics/tensors
All other classes will be used and updated in future or removed.

### Next steps ###
Adding the homogeneouse reactor calculator and all relevant classes

### How to get the repository and work with it ###
* Feel free to compile it where ever you want, but normally its nice to have a fixed folder for _user compiled stuff_
* Make a new folder

```bash
mkdir -p $HOME/yourLocation
```

* Switch to the new folder

```bash
cd $HOME/yourLocation
```

* Clone the repository to the new folder

```bash
git clone https://github.com/shor-ty/thermoKineticCalculator.git thermoKineticCalculator
```

* Switch to the repository directory

```bash
cd thermoKineticCalculator
```

* Load the environment variables

```bash
source environment/bash
```

* The environment should be loaded everytime you want to use the applications. Thus, an alias is recommended or your simply load it each time a terminal is opened.

```bash
echo "alias tkc="source $HOME/yourLocation/thermoKineticCalculator/environment/bashrc" >> $HOME/.bash_aliases
```

* Compile the libraries

```bash
cd $TKC_PROJECT_DIR/src
make
```

* Compile the applications

```bash
cd $TKC_PROJECT_DIR/application/homogeneousReactor
make
```

* You are ready to start with the tutorials

```bash
cd $TKC_PROJECT_DIR/tutorials
```
