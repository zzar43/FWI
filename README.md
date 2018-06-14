# FWI

This repo contains some notes and codes about full waveform inversion.

Here is a note about [adjoint equation](https://zzar43.github.io/FWI/notes/adjoint_eq.html) based on elastic wave equation and least square objective function.

## Idea
- Using struct to describe velocity model, source and receiver information
- Wavefield data and recorded data should be designed for multi-frequency and multi-sources directly
- Input and output format should be vector. Corresponding to the discrete math form
- Function should be treated as functional operator

## Unit

- Velocity: m/s
- Density: kg/m^3
- Sampling frequency: Fs (Hz)
- Time: t (s)
