# Interactive ARAP

***As-rigid-as-possible surface deformation* (ARAP)** project for the 3D Scanning and Motion Capture course WS2022-2023. This project is a re-implementation of Sorkine et al [[_1_]](https://igl.ethz.ch/projects/ARAP/arap_web.pdf). Our motivation for choosing this topic was that it is interactive and independent of additional hardware. In the [***report***](/Report.pdf), we review related work, the methods used for the algorithm, our results, and the conclusion, which includes the challenges we encountered, and our future work.

## SETUP

**Requirements**
- OpenGL
- all other requirements are installed via the Cmake file
- *(Optional)* You can run the project parallel with OpenMP

Use CMAKE to generate Visual Studio solution with the following flags

![Config Flags](/image.png)

## RUN

1. Open the generated solution in Visual Studio
2. Set Startup Project to interactive_arap in the Solution Explorer

Note: If you get the error "MSB3073	The command setlocal..." while building, go to interactive_arap Project Properties -> Build Events -> Post Build Events -> Delete set local command, apply and build again.

### Contributors

- [Andrea Solanas de Vicente](mailto:andrea.solanasvicente@gmail.com)
- [Ankur Deria](mailto:ankurderia1999@gmail.com)
- [Bendeguz Timar](mailto:timar.bendi@gmail.com)
- [Michael Dey](mailto:micha.eldey@yahoo.com)

### References

[[_1_]](https://igl.ethz.ch/projects/ARAP/arap_web.pdf) Olga Sorkine and Marc Alexa. As-rigid-as-possible surface modeling. In *Proceedings of the Fifth Eurographics Sym- posium on Geometry Processing*, SGP ’07, page 109–116, Goslar, DEU, 2007. Eurographics Association.
