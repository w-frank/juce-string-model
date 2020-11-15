# JUCE String Model

Realistic string instrument model using physical modelling. Includes a delay line
to create echo patterns in the stereo sound field. Built using the [JUCE](https://github.com/juce-framework/JUCE) 
audio application framework.

Based on the excellent tutorial by [JUCE](https://docs.juce.com/master/tutorial_dsp_delay_line.html).

## Build Environment
* JUCE v6.0.4
* Ubuntu 18.04.5 LTS

## Build
```bash
cd Builds/LinuxMakefile
make
```

## Run
```bash
cd Builds/LinuxMakefile/build
./StringModel
```