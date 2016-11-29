# BrownDyeTools

[BrownDye](http://browndye.ucsd.edu) plugin for [Pymol](http://ww.pymol.org)

This Pymol plugin allows setting up, running and analyzing Brownian
dynamics simulation using the BrownDye collection of applications.

![main window](https://github.com/rokdev/BrownDyeTools/blob/master/main.png)

## Requirements

* [PDB2PQR](http://www.poissonboltzmann.org) Download and install the
  precompiled binary version of PDB2PQR.

* [APBS](http://www.poissonboltzmann.org) Download and install the
  precompiled binary version of APBS.

* [BrownDye](http://browndye.ucsd.edu) There is no precompiled binary
  distribution, you have to download and compile the source code.

Currently only Linux OS is supported. MacOS is untested.

## Installation

* Download and install [pymol](http://sourceforge.net/projects/pymol/)

* Download [BrownDyeTools](https://github.com/rokdev/BrownDyeTools)

* Open `PyMol` and install `BrownDyeTools`: `Plugin -> Manage Plugins -> Install` then 
  locate `BrownDyeTools.py` file.

## Using

* Start `pymol` and open BrownDyeTools plugin: `Plugin -> BrownDye Tools`

### Configuration Tab

* Create a randomly named project directory or browse for a custom
  directory.

* Set `PDB2PQR_PATH`, `APBS_PATH` and `BD_PATH` locations. If you have these set as environment variables in your shell they will be picked up by the plugin.

### PQR files Tab

* Pick molecule 0 and 1 PDB files, or use Pymol selection.

* Press "Create PQR files" to create PQR files from the PDB files using pre-set force field.

* Or load your own PQR files.

### APBS Tab

* Press `Calculate grid size` for both molecule 0 and 1 to get
  automatic grid dimensions or fill in your customized grid size
  values.

* Adjust APBS options if necessary and press `Run APBS to generate
  grids` to launch two APBS calculation (for molecule 0 and 1). This
  might take a couple of minutes, depending on the molecule size and
  your computing hardware.

### Reaction criteria Tab

* You can use either the default protein-protein contact file or load
  your own.

* Adjust the reaction criteria if desired and press `Create
  reaction files` to generate the necessary reaction input files.
  

### BD setup Tab

* Adjust - if necessary - BrownDye input file parameters. For detailed
  description of the input parameters see BrownDye documentation. 
  `Get Debye length` will automatically fetch the APBS calculated value.

* Press `Generate BD input files` to create the necessary BrownDye
  input files. This might take a couple of minutes, depending on the
  size of molecules 0 and 1.

### BD simulation Tab

* You can choose either to run BrownDye simulation within the plugin
  or send it to background. Obviously, the foreground (in the plugin)
  execution is appropriate for shorter simulation. When running in the
  plugin various statistics are updated during the run and displayed.

### Analysis Tab

* After simulation is finished you can pick and analyze the created
  trajectory files. The `Analyze` button will print how many
  successful events were recorder in the particular trajectory files
  and their index numbers. You can then pick one successful trajectory
  (`Select trajectory file`) and convert it to XYZ file (`Convert to
  xyz trajectory`) and visualize it in Pymol (`Load xyz trajectory`).


****

## License

This software is released under the terms of GNU GPL3 license.
For more details please see the accompanying documentation.

(c) 2016 [National Biomedical Computation Resource](http://nbcr.ucsd.edu)
