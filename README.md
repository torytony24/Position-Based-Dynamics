# Position Based Dynamics

This repository provides reference implementations of
- **Position Based Dynamics (PBD)**

    _Müller, Matthias, et al. "Position based dynamics." Journal of Visual Communication and Image Representation 18.2 (2007): 109-118._

- **Extended Position Based Dynamics (XPBD)**

    _Macklin, Miles, et al. "XPBD: Position-based simulation of compliant constrained dynamics." Proceedings of the 9th ACM SIGGRAPH/Eurographics Symposium on Computer Animation. 2016._

The project is done with the help of [Graphics and Media Lab](http://graphics.snu.ac.kr/) (Department of Electrical and Computer Engineering, Seoul National University).

## ✨ Features

### PBD implementation

![Image](https://github.com/user-attachments/assets/406e113c-ff4e-494d-aad0-b13a7eee7266)

### XPBD implementation

![Image](https://github.com/user-attachments/assets/2f44e207-767e-4a2f-8973-336ce9ca3b6e6)

### Object example details

- PBD_v1 / XPBD
    - Cloth simulation with `clothMesh.obj`
    - Can change the fixed vertices
- PBD_v2
  - Volumetric object simulation with `50Cube.obj` and `35Sphere.obj`

## 🎥 Demo

HERE

## 🔍 Technical Details

For further informations, refer to the original papers or the presentation files [PBD_review_240118.pdf](./PBD_review_240118.pdf) and [XPBD_review_240207.pdf](./XPBD_review_240207.pdf).

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
