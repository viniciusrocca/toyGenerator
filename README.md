# toyGenerator
A simplified Monte carlo event generator used to simulate the generation of a Higgs boson as a signal and your decay in two photons.

## Generator Structure 

The following image represents the structure of the simplified Monte Carlo simulator. Where:

  * Particles 1 and 2 are partons.
  * Particle X is a Higgs boson.
  * Particle Y is a parton.
  * Particles A and B are photons.
  * The squares represents the detectors.

![Image of Yaktocat](https://github.com/viniciusrocca/toyGenerator/blob/master/generatorStructure.jpg)


The generator is divided in four blocks:

  * Initial state block: Sampling of the partons energys according to a PDF.
  * Collision block: Simulate the collision of the partons generating X (Higgs) and Y (parton).
  * Decay block: Simulate the decay of X (Higgs) in A and B (two photons).
  * Detector block: Simulate the smearing effect caused by the detector.

## Running

A basic example to run this generator is:

```
./EventGenerator.py -n 1000
```

To run 1000 events of Higgs boson generation as a singal and his decay in two photons.


## Possible inputs

 * Number of events ```-n <n_events>```.
 * Proton energy at the center of mass in TeV ```-e <proton_energy>```.
 * Mass of particle X in GeV ```-x <x_mass>```
 * Mass of particle Y in MeV ```-y <y_mass>```
 * Mass of parton 1 in MeV ```-p <parton1_mass>```
 * Mass of parton 2 in Mev ```-q <parton2_mass>```
 * Mass of particle A in GeV ```-a <a_mass>```
 * Mass of particle B in GeV ```-b <b_mass>```


