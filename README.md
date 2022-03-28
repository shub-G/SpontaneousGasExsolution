# SpontaneousGasEvolution

This repository contains source code for the g-cubed manuscript titled "Spontaneously exsolved free gas during major storms as an ephemeral gas source for pockmark formation".

## Plain language summary of the study
Thousands of pockmarks, circular morphological depressions on the seafloor, were reported in south-eastern North Sea, presumably formed in response to wave motions during major storms.

It has been hypothesized that the pockmarks formed as pre-existing shallow free gas pockets were mobilized by pressure changes of the waves. However, the mechanisms that could have mobilized free gas are not yet constrained. Moreover, large scale free gas accumulations have not been reported in this region, and therefore, commonly invoked mechanisms like tensile failure and breaching of capillary seals are hard to justify as they rely on the presence of pre-existing gas pockets. 

Here, through modelling studies, we tackle the question of the source of the observed free gas. Our study consists of two parts: 
First, assuming that some hitherto unknown shallow free-gas pocket is indeed present, we test whether storm-induced pressure changes could be sufficient to breach capillary seals. We find that free gas damps the pressure changes due to its high compressibility, making the mobilization of pre-existing gas unlikely. 
In the second part, we propose an alternative mechanism where free gas spontaneously appears due to exsolution from pore-fluids, primarily driven by the pressure-sensitivity of gas solubility.
We test the feasibility of this mechanism using an idealized test setting, and show how periodic pressure changes can result in the appearance of a persistent gas phase that could offer a possible explanation for the elusive gas source linked to these pockmarks.

## Numerical settings

1. 1D scenarios are simulated with a pre-existing free-gas pocket burried below an idealized capillary barrier. These simulations test whether gas can be mobilized past a capillary seal purely through pumping by storm-related pressure changes, as proposed in the original work of Krämer et al. (2017).

2. 1D scenarios are simulated *without* any pre-existing free gas pockets. These simulations test the main hypothesis of *spontaneous appearance of free gas by exsolution from pore-fluids*, driven by storm-related pressure changes.

3. 2D scenarios are simulated to show the impact of lateral gradients on this proposed mechanism.

4. An example 2D scenario is simulated with evolving permeability as a function of changes in pore-pressure. This shows the emergence of periodic structures in the sediment.

# How to use this source code?

This code requires DUNE-PDELab version 2.8, which is an open-source toolbox for solving systems of partial differential equations, available freely at: https://gitlab.dune-project.org/pdelab/

