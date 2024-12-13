# Volcanic Skies
This is the code used for generating the results shown in [Volcanic Skies](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.15034) \[Pretorius, P. C., Gain, J., Lastic, M., Cordonnier, G., Chen, J., Rohmer, D., & Cani, M. P. (2024, May). Volcanic Skies: coupling explosive eruptions with atmospheric simulation to create consistent skyscapes. In Computer Graphics Forum (Vol. 43, No. 2, p. e15034).\]

It requires the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [OpenVDB](https://www.openvdb.org/) libraries to be installed. The Eigen header only files can be placed in `lib\Eigen\` for easy compilation.

The MergeTestbed executable (`src\Tests\MergeTest.cpp`) contains the scenarios described in the paper. In particular, the capWilson, skirts3, and AshRainBigFinal scenarios are described in the paper.
