# magnetictb-abstract-group

This crate implements the frozen table-based finite-group API from
`magnetictb.abstract_group.table_fixture` version 1.

The public boundary uses zero-based ordered indices:

- `multiplication_table[left][right] = left * right`;
- `action_table[operation][source] = target`;
- right cosets are `representative * subgroup`;
- left cosets are `subgroup * representative`;
- Schreier decomposition satisfies
  `operation * representative[source] = representative[target] * h`.

`GroupEnumeration` is intentionally not implemented. Its opaque
callback-driven elements are explicitly excluded from fixture schema version
1.

The crate has no Mathematica runtime dependency. Mathematica is used only to
produce the immutable local cross-language fixtures.
