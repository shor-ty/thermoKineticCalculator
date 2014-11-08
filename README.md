# README #

This README  includes the steps which are necessary to get this application up and running.

### OpenFOAM Versions ###
* Generated for all versions (not a OpenFOAM utility)

### What is this repository for? ###
* This application calculate the laminar flamelets for your individual chemistry
* You can use either the laminar flamelets or extend this to the turbulent flamelets with distribution functions
* It is the pre-processing step when using the flamelet model that I update (build by Alberto Cuoci)
* Version 1.0
* Developed by [Holzmann-cfd](https://holzmann-cfd.de)

### Prerequisists ###
* The following utilitys are necessary for successful compiling

### How do I get set up? ###
* Feel free to compile it where ever you want, but normally its nice to have a fixed folder for _user compiled stuff_
* Make a new folder
> mkdir -p $FOAM_RUN/../OpenFOAM_extensions
* Switch to the new folder
> cd $FOAM_RUN/../OpenFOAM_extensions
* Clone the repository to the new folder
> git clone https://shor-ty@bitbucket.org/shor-ty/flameletcreator.git
* Switch to the repository directory
> cd flameletcreator
* Compile the application
> g++
* Finished

### Contribution guidelines ###
* If you have questions, hints or any suggestions please email me to Tobias.Holzmann@Holzmann-cfd.de

### Other stuff ###
* Thanks to Oliver Borm for some hints