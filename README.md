# SURFACE-Bind
Surfaceome Targetability for Functional Binder Design

## The workflow
We computationally mined the entire human surfaceome, predicting targetable binding interfaces for protein binder design strategies, via assigning a set of geometric, and chemical scores. These unbound-state scoring led to an average 5 target sites per protein for over 2,800 protein entries in the surfaceome dataset. To further evaluate the targetability of these binding interfaces, we utilized a library of 640,000 continuous structural fragments (seeds) with different secondary structure elements and performed over 3 billion pairwise docking, scoring, sorting, and selection which provided a bound-state score of the target sites, as well as a list of high-quality candidate seeds. In result, we are reporting the most complete set of targetable binding interfaces on the human surfaceome, linked to high-quality seeds that can initiate further protein-based therapeutic designs.

## Authors and acknowledgment
This project is developed through collaborations between EPFL, Novo Nordisk, and Inria.

## License
Please see the License file.

## Citations
@article{balbi2024mapping,
  title={Mapping targetable sites on the human surfaceome for the design of novel binders},
  author={Balbi, Petra EM and Sadek, Ahmed and Marchand, Anthony and Yu, Ta-Yi and Damjanovic, Jovan and Georgeon, Sandrine and Schmidt, Joseph and Fulle, Simone and Yang, Che and Khakzad, Hamed and others},
  journal={bioRxiv},
  pages={2024--12},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}