# INTERACTIVE ARAP

***As-rigid-as-possible surface deformation* (ARAP)** project for the 3D Scanning and Motion Capture course WS2022-2023. This project is a re-implementation of Sorkine et al [[_1_]](https://igl.ethz.ch/projects/ARAP/arap_web.pdf). Our motivation for choosing this topic was that it is interactive and independent of additional hardware. In the [***report***](/Report.pdf), we review related work, the methods used for the algorithm, our results, and the conclusion, which includes the challenges we encountered, and our future work.

### LOCAL ENVIRONMENT SETUP

You can easily download and run the project without installing a lot of libraries. Please use the following instructions to do so.

- You need to have **OpenGL** installed on your system to build and run the project.

- The library, **[libigl](https://libigl.github.io/tutorial/)**, will be installed by itself when using the given CmakeLists file.

- Use CMAKE to generate Visual Studio solution. Use the following flags in cmake config ![Config Flags](/image.png)

- *(Optional)* If you have **OpenMP** installed on your system, the project will run in parallel.

### HOW TO RUN

- Build and run the generated solution in Visual Studio after setting interactive_arap as Startup Project in the Solution Explorer. (Use Release version for better performance)

Note: If you get the error "MSB3073	The command setlocal..." while building, go to interactive_arap Project Properties -> Build Events -> Post Build Events -> Delete set local command and apply  and build again.

### Contributors

- Andrea Solanas de Vicente
- Michael Dey
- Ankur Deria
- Bendeguz Timar

### References

[[_1_]](https://igl.ethz.ch/projects/ARAP/arap_web.pdf) Olga Sorkine and Marc Alexa. As-rigid-as-possible surface modeling. In *Proceedings of the Fifth Eurographics Sym- posium on Geometry Processing*, SGP ’07, page 109–116, Goslar, DEU, 2007. Eurographics Association.
