# sla-planet-sim
Simulating the birth of planets using semi-Lagrangian advection and spectral methods.

Based on the research and code contained in 
SIMULATING THE BIRTH OF PLANETS: A SPECTRAL SEMI-LAGRANGIAN HYDRODYNAMIC APPROACH,
a Master's thesis by Wendy Crumrine. Code ported to C from MATLAB for parallelization.

[Introduction](#introduction) •
[Getting started](#getting-started) •
[Installation](#installation) •
[Configuration](#configuration) •
[TODOs](#todos) •
[Integrations](#third-party-integrations)

## Introduction
Planetary formation takes place in a disc of dust and gas rotating around a star.  How planets grow from meter-sized boulders to kilometer-sized planetessimals is still a mystery, as our current models suggest that drag forces within the gas should cause these relatively small bodies to spiral into the star.  

In this project we address the so-called "meter-sized gap problem" by modeling gas and dust as fluids, using a combination of Eulerian and Lagrangian dynamics.  A key component to our approach is quantifying the vorticity of the gas and dust at every point in our sample space.  The underlying model includes five PDE's which are solved using spectral methods.  

## Installation

|![](images/time-lapse.gif)|
|:--:| 
| Simulation of changing dust densities in proto-planetary disc leading to the formation of planetesimals. |

## Getting started

jasldfkjdsjfk

## TODOs
|![](images/time-lapse.gif)|
|:--:| 
| Simulation of changing dust densities in proto-planetary disc leading to the formation of planetesimals. |

fft.h / fft.c - Fast fourier transform functions based on FFTW3. 
sla.h / sla.c - Variable declarations and function defintions. 
