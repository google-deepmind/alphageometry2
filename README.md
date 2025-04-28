# Gold-medalist Performance in Solving Olympiad Geometry with AlphaGeometry2

We present AlphaGeometry2, a significantly improved version of AlphaGeometry
introduced in [Trinh et al. (2024)](https://www.nature.com/articles/s41586-023-06747-5),
which has now surpassed an average gold medalist in solving Olympiad geometry
problems.

This repository contains code necessary to reproduce solving some olympiad
geometry problems (for example, IMO 2005 P1, IMO 2008 P6, IMO 2013 P3 and
others).
This is done by running the symbolic part of AlphaGeometry2 system called DDAR.
DDAR runs on problems formalized into the AlphaGeometry2 language and
supplemented by the point coordinates. Some problems can be solved by DDAR
alone and others require auxiliary points provided by a language model.

## Installation

Installation is done in a virtual environment:

```
$ sudo apt install python3
$ sudo apt install python3-pip
$ sudo apt install python3-venv
$ python3 -m venv ag2
$ source ag2/bin/activate
$ pip install numpy
```

## Usage

Command below runs the prover on some IMO problems solved by AlphaGeometry2:

```
$ python -m test
```

You will see the following output:

```

We run the logical core DDAR on some challenging IMO problems, with manually
provided auxiliary points. This is not a full AlphaGeometry system, only a
test of the logical core.

Problem: 2001_p5a
......... Proven :-)
Problem: 2003_p4b
..... Proven :-)
Problem: 2005_p1
...... Proven :-)
Problem: 2008_p1b
........ Proven :-)
Problem: 2008_p6
...... Proven :-)
Problem: 2009_p4a
..... Proven :-)
Problem: 2009_p4b
...... Proven :-)
Problem: 2011_p6
...... Proven :-)
Problem: 2013_p3
....... Proven :-)
Problem: 2018_p6
..... Proven :-)
Problem: 2019_p2
.... Proven :-)
Problem: 2021_p3
......... Proven :-)
Problem: 2023_p6
....... Proven :-)
```

## Citing this work

Please cite this work as:

```
@misc{chervonyi2025goldmedalistperformancesolvingolympiad,
      title={Gold-medalist Performance in Solving Olympiad Geometry with AlphaGeometry2},
      author={Yuri Chervonyi and Trieu H. Trinh and Miroslav Olšák and Xiaomeng Yang and Hoang Nguyen and Marcelo Menegali and Junehyuk Jung and Vikas Verma and Quoc V. Le and Thang Luong},
      year={2025},
      eprint={2502.03544},
      archivePrefix={arXiv},
      primaryClass={cs.AI},
      url={https://arxiv.org/abs/2502.03544},
}
```

## License and disclaimer

Copyright 2025 Google LLC

All software is licensed under the Apache License, Version 2.0 (Apache 2.0);
you may not use this file except in compliance with the Apache 2.0 license.
You may obtain a copy of the Apache 2.0 license at:
https://www.apache.org/licenses/LICENSE-2.0

All other materials are licensed under the Creative Commons Attribution 4.0
International License (CC-BY). You may obtain a copy of the CC-BY license at:
https://creativecommons.org/licenses/by/4.0/legalcode

Unless required by applicable law or agreed to in writing, all software and
materials distributed here under the Apache 2.0 or CC-BY licenses are
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
either express or implied. See the licenses for the specific language governing
permissions and limitations under those licenses.

This is not an official Google product.
