# Position Based Dynamics

This repository provides reference implementations of
- **Position Based Dynamics (PBD)**

    _Müller, Matthias, et al. "Position based dynamics." Journal of Visual Communication and Image Representation 18.2 (2007): 109-118._

- **Extended Position Based Dynamics (XPBD)**

    _Macklin, Miles, et al. "XPBD: Position-based simulation of compliant constrained dynamics." Proceedings of the 9th ACM SIGGRAPH/Eurographics Symposium on Computer Animation. 2016._

The project is done with the help of [Graphics and Media Lab](http://graphics.snu.ac.kr/) (Department of Electrical and Computer Engineering, Seoul National University).

## ✨ Features

### PBD implementation

<img width="1280" height="720" alt="PBD" src="https://github.com/user-attachments/assets/32a57f58-5f8d-424e-be43-cf1d7aabfeda" />

### XPBD implementation

<img width="1280" height="720" alt="XPBD" src="https://github.com/user-attachments/assets/91b9fb59-965e-49b5-bd4d-1e7ae3c81738" />

### Object example details

- PBD_v1 / XPBD
    - Cloth simulation with `clothMesh.obj`
    - Can change the fixed vertices
- PBD_v2
  - Volumetric object simulation with `50Cube.obj` and `35Sphere.obj`

## 🎥 Demo

https://github.com/user-attachments/assets/adac0578-1fb2-4b76-8c8d-6888915881c6

## 🔍 Technical Details

For further informations, refer to the original papers or the presentation files [PBD_review_240118.pdf](./PBD_review_240118.pdf) and [XPBD_review_240207.pdf](./XPBD_review_240207.pdf).

Research and study archive for the project can be found on [this site](https://researchandstudy.notion.site/Cloth-Simulation-cb04296d17d346e0b63926a6b159f2e6?pvs=143).

## ▶️ How to Run

### Requirements
- C++17 or later
- CMake (≥ 3.10)

### Run PBD
```bash
# Build using CMake
git clone https://github.com/torytony24/Position-Based-Dynamics.git
cd Position-Based-Dynamics/PBD_v2/PBD/PBD
mkdir build && cd build
cmake ..
make

# Run main.cpp
./PBD
```

### Run XPBD
```bash
# Build using CMake
git clone https://github.com/torytony24/Position-Based-Dynamics.git
cd Position-Based-Dynamics/XPBD/PBD/PBD
mkdir build && cd build
cmake ..
make

# Run main.cpp
./PBD
```

## 📝 References

- [_Müller, M., Heidelberger, B., Hennix, M., & Ratcliff, J. (2007). Position Based Dynamics._](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf)
- [_Macklin, M., Müller, M., Chentanez, N., & Jeschke, S. (2016). XPBD: Position-based simulation of compliant constrained dynamics._](https://matthias-research.github.io/pages/publications/XPBD.pdf)
