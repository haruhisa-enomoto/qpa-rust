# Pinned QPA manual chapters

The original QPA manual remains the most detailed mathematical documentation
for many specialized operations. The files below are pinned with the source
used for this port. They use GAP syntax; combine them with the generated
[QPA-to-Python operation index](../python/qpa-operation-index.md) or the
[GAP/QPA-to-Rust map](gap-to-rust.md).

| Topic | Pinned QPA manual source |
| --- | --- |
| Introduction and aims | [`chapter_introduction.xml`](../../upstream/qpa/doc/chapter_introduction.xml) |
| Original QPA quickstart | [`chapter_quickstart.xml`](../../upstream/qpa/doc/chapter_quickstart.xml) |
| Quivers and paths | [`chapter_quivers.xml`](../../upstream/qpa/doc/chapter_quivers.xml) |
| Path algebras and quotients | [`chapter_path_algebras.xml`](../../upstream/qpa/doc/chapter_path_algebras.xml) |
| Right modules | [`chapter_right_modules.xml`](../../upstream/qpa/doc/chapter_right_modules.xml) |
| Module homomorphisms | [`chapter_homomorphisms.xml`](../../upstream/qpa/doc/chapter_homomorphisms.xml) |
| Homological algebra | [`chapter_homological_algebra.xml`](../../upstream/qpa/doc/chapter_homological_algebra.xml) |
| Auslander--Reiten theory | [`chapter_AR-theory.xml`](../../upstream/qpa/doc/chapter_AR-theory.xml) |
| Chain complexes | [`chapter_chain_complexes.xml`](../../upstream/qpa/doc/chapter_chain_complexes.xml) |
| Projective and derived computations | [`chapter_proj_and_derived.xml`](../../upstream/qpa/doc/chapter_proj_and_derived.xml) |
| Gröbner bases | [`chapter_groebner_bases.xml`](../../upstream/qpa/doc/chapter_groebner_bases.xml) |
| QPA's GBNP interface | [`chapter_gbnp.xml`](../../upstream/qpa/doc/chapter_gbnp.xml) |
| Combinatorial representations | [`chapter_combinatorialrep.xml`](../../upstream/qpa/doc/chapter_combinatorialrep.xml) |
| Degeneration order and finite type | [`chapter_degorderfinitetype.xml`](../../upstream/qpa/doc/chapter_degorderfinitetype.xml) |

These chapters describe QPA's intended operations. When the pinned
implementation differs, the executable source, oracle fixtures, and recorded
[deviations](../../porting/deviations.md) are authoritative for this port.
