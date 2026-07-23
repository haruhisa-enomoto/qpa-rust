# Notices and attribution

This repository is a modified, source-faithful Rust translation of QPA and of
the part of GBNP transitively used by QPA. It is not an official release of
QPA, GBNP, GAP, or Numerica. Copyright in contributed Rust code remains with
its contributors; copyright in translated source material remains with the
respective upstream authors.

The workspace is distributed under `GPL-2.0-or-later`. See `LICENSE`.

## QPA

QPA (Quivers and Path Algebras) is the primary computational and API source for
this port.

- Project: <https://github.com/gap-packages/qpa>
- Pinned version/revision: 1.37 / `8d5f847258d4853700ce3a5d947db2a1441109f1`
- License: GPL-2.0-or-later
- Upstream license and author information: `upstream/qpa/LICENSE` and
  `upstream/qpa/PackageInfo.g` in a development checkout

## GBNP

The Rust GBNP crate translates the transitive noncommutative Groebner machinery
that QPA invokes.

- Project: <https://github.com/gap-packages/gbnp>
- Pinned version/revision: 1.1.0 / `e0b46c1aac16f3fbabe743e3169dc6bd86bdc654`
- License: LGPL-2.1-or-later
- Principal upstream copyright notice: Copyright 2001-2010 Arjeh M. Cohen,
  Die A.H. Gijsbers, Jan Willem Knopper, and Chris Krook

## GAP

GAP supplies dispatch, field, matrix, basis, and generic algebra semantics used
as a transitive specification and by the differential oracle. The pinned GAP
source and binary are development-only references and are not linked into the
Rust crates. The compact `qpa-core` Conway table is generated from GAP's
GPL-2.0-or-later `lib/conwdat1.g`, `conwdat2.g`, and `conwdat3.g`; its generator
records the pinned source revision in the generated file.

- Project: <https://github.com/gap-system/gap>
- Pinned version/revision: 4.17dev / `4efe877859ab0bc0ac041bf635a996fb897ee0c9`
- License: GPL-2.0-or-later

## Numerica

Numerica provides exact scalar and linear-algebra primitives. It does not
replace QPA algorithms.

- Project: <https://github.com/symbolica-dev/numerica>
- Pinned crate version: 2.1.0
- Inspected revision: `da17228d1d2a154da06efbfcda1c739770192f25`
- License: MIT

The machine-readable authoritative dependency and revision record is
`porting/upstream-lock.toml`. Development clones in `upstream/` are excluded
from source packages and future Python wheels.
