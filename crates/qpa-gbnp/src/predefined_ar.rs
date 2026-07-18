//! QPA's predefined finite-type numerical Auslander--Reiten quivers.
//!
//! The tables and parameter formulas translate QPA 1.37
//! `lib/degorderfinitetype.gi`, `PredefARQuivers`, lines 1303--1794.

use thiserror::Error;

use crate::{
    ArMesh, ArOrbit, ArOrbitPosition, ArQuiverDescription, ArQuiverNumerical, DegenerationError,
};

/// The predefined numerical AR-quiver selectors accepted by QPA.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum PredefinedArQuiver {
    /// The subspace orientation of Dynkin type D4.
    D4Subspace,
    /// The subspace orientation of Dynkin type E6.
    E6Subspace,
    /// Equioriented A3 modulo its length-two path.
    A3Zero,
    /// QPA's `A1,2 zero` example.
    A12Zero,
    /// QPA's `A2,2 comm` example.
    A22Commutative,
    /// QPA's `R nilp` example.
    RNilpotent,
    /// An entry from the Bongartz--Gabriel list, with an optional parameter.
    BongartzGabriel {
        /// List number.
        number: usize,
        /// Required for list numbers 1, 3, 4, 6, and 10; absent otherwise.
        parameter: Option<usize>,
    },
    /// One of QPA's alternative compact orbit encodings (`BGo`).
    BongartzGabrielOrbits {
        /// List number; QPA supplies 2, 13, and 14.
        number: usize,
    },
}

/// Raw data returned by QPA's `PredefARQuivers` before construction.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PredefinedArQuiverData {
    /// Number of indecomposable modules.
    pub number_of_indecomposables: usize,
    /// Number of indecomposable projective modules.
    pub number_of_projectives: usize,
    /// Explicit-mesh or compact-orbit description.
    pub description: ArQuiverDescription,
}

impl PredefinedArQuiverData {
    /// Apply QPA's `ARQuiverNumerical` constructor to this predefined data.
    ///
    /// # Errors
    ///
    /// Propagates numerical AR-matrix validation failures.
    pub fn into_ar_quiver(self) -> Result<ArQuiverNumerical, DegenerationError> {
        ArQuiverNumerical::new(
            self.number_of_indecomposables,
            self.number_of_projectives,
            self.description,
        )
    }
}

/// Failure of QPA's predefined-name/parameter dispatch.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PredefinedArQuiverError {
    /// QPA has no branch for the requested selector.
    #[error("wrong number or type of arguments")]
    Unsupported,
}

/// QPA's stable output for `PredefARQuivers("what")`.
pub const PREDEFINED_AR_QUIVER_CATALOG: &str = concat!(
    "1) (\"BG\",i), for i=2,5,7,8,9,11,12,13,14, algebra no. i from Bongartz-Gabriel list ",
    "of maximal finite type with 2 simples;\n",
    "2) (\"BG\",i,r), for i=1,3,4,6,10; r>=1; algebra no. i with parameter r ",
    "from Bongartz-Gabriel list of maximal finite type with 2 simples;\n",
    "3) (\"D4 subspace\"), path algebra of Dynkin quiver D4 with subspace\n orientation of arrows;\n",
    "4) (\"E6 subspace\"), path algebra of Dynkin quiver E6 with subspace\n orientation of arrows;\n",
    "5) (\"A3 zero\"), path algebra of equioriented Dynkin quiver A3  modulo unique zero relation of length 2;\n",
    "6) (\"A1,2 zero\"), path algebra of Euclidean quiver A1,2 modulo zero relation  ",
    "(i.e. Quiver(3,[[3,1,\"a\"],[3,2,\"b\"],[2,1,\"c\"]])/(bc) );\n",
    "7) (\"A2,2 comm\"), path algebra of Euclidean quiver A2,2 modulo unique commutativity relation;\n",
    "8) (\"R nilp\"), one arrow + one loop modulo nilpotency deg. 2  ",
    "(i.e. Quiver(2,[[2,1,\"a\"],[1,1,\"b\"]])/(b^2) );\n",
);

/// Translate QPA's complete `PredefARQuivers` data dispatch.
///
/// # Errors
///
/// Returns QPA's wrong-arguments failure for an unsupported number, missing or
/// superfluous parameter, or zero parameter outside the documented `r >= 1`
/// domain.
pub fn predefined_ar_quiver_data(
    selector: PredefinedArQuiver,
) -> Result<PredefinedArQuiverData, PredefinedArQuiverError> {
    match selector {
        PredefinedArQuiver::D4Subspace => Ok(d4_subspace()),
        PredefinedArQuiver::E6Subspace => Ok(e6_subspace()),
        PredefinedArQuiver::A3Zero => Ok(a3_zero()),
        PredefinedArQuiver::A12Zero => Ok(a12_zero()),
        PredefinedArQuiver::A22Commutative => Ok(a22_commutative()),
        PredefinedArQuiver::RNilpotent => Ok(r_nilpotent()),
        PredefinedArQuiver::BongartzGabriel {
            number,
            parameter: Some(parameter),
        } if parameter > 0 => bg_parameterized(number, parameter),
        PredefinedArQuiver::BongartzGabriel {
            number,
            parameter: None,
        } => bg_fixed(number),
        PredefinedArQuiver::BongartzGabrielOrbits { number } => bg_orbits(number),
        PredefinedArQuiver::BongartzGabriel { .. } => Err(PredefinedArQuiverError::Unsupported),
    }
}

/// Construct one of QPA's complete predefined numerical AR quivers.
///
/// # Errors
///
/// Propagates selector and numerical-construction failures.
pub fn predefined_ar_quiver(
    selector: PredefinedArQuiver,
) -> Result<ArQuiverNumerical, PredefinedArQuiverConstructionError> {
    Ok(predefined_ar_quiver_data(selector)?.into_ar_quiver()?)
}

/// Failure while selecting or constructing a predefined numerical AR quiver.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PredefinedArQuiverConstructionError {
    /// The predefined selector is unsupported.
    #[error(transparent)]
    Selector(#[from] PredefinedArQuiverError),
    /// QPA's numerical constructor rejected its predefined data.
    #[error(transparent)]
    Construction(#[from] DegenerationError),
}

macro_rules! rows {
    ($( [$($entry:expr),* $(,)?] ),* $(,)?) => {
        vec![$(vec![$($entry),*]),*]
    };
}

fn mesh_data(projectives: usize, rows: Vec<Vec<usize>>) -> PredefinedArQuiverData {
    let number_of_indecomposables = rows.len();
    let description = rows
        .into_iter()
        .map(|mut row| {
            let tau = row.pop().expect("QPA mesh rows are nonempty");
            ArMesh::new(
                row.into_iter().map(|position| position - 1).collect(),
                (tau != 0).then_some(tau.saturating_sub(1)),
            )
        })
        .collect();
    PredefinedArQuiverData {
        number_of_indecomposables,
        number_of_projectives: projectives,
        description: ArQuiverDescription::Meshes(description),
    }
}

fn orbit(length_minus_one: usize, predecessors: &[(usize, usize)]) -> ArOrbit {
    ArOrbit::new(
        length_minus_one,
        predecessors
            .iter()
            .map(|&(orbit, shift)| ArOrbitPosition::new(orbit - 1, shift))
            .collect(),
    )
}

fn orbit_data(projectives: usize, orbits: Vec<ArOrbit>) -> PredefinedArQuiverData {
    let number_of_indecomposables = orbits.iter().map(|orbit| orbit.length_minus_one + 1).sum();
    PredefinedArQuiverData {
        number_of_indecomposables,
        number_of_projectives: projectives,
        description: ArQuiverDescription::Orbits(orbits),
    }
}

fn bg_parameterized(
    number: usize,
    parameter: usize,
) -> Result<PredefinedArQuiverData, PredefinedArQuiverError> {
    match number {
        1 => Ok(bg1(parameter)),
        3 => Ok(bg3(parameter)),
        4 => Ok(bg4(parameter)),
        6 => Ok(bg6(parameter)),
        10 => Ok(bg10(parameter)),
        _ => Err(PredefinedArQuiverError::Unsupported),
    }
}

fn bg_fixed(number: usize) -> Result<PredefinedArQuiverData, PredefinedArQuiverError> {
    match number {
        2 => Ok(bg2()),
        5 => Ok(bg5()),
        7 => Ok(bg7()),
        8 => Ok(bg8()),
        9 => Ok(bg9()),
        11 => Ok(bg11()),
        12 => Ok(bg12()),
        13 => Ok(bg13()),
        14 => Ok(bg14()),
        _ => Err(PredefinedArQuiverError::Unsupported),
    }
}

fn bg_orbits(number: usize) -> Result<PredefinedArQuiverData, PredefinedArQuiverError> {
    let orbits = match number {
        2 => vec![
            orbit(7, &[(2, 2)]),
            orbit(9, &[(3, 0)]),
            orbit(15, &[(4, 15)]),
            orbit(15, &[(3, 0), (5, 7)]),
            orbit(7, &[(4, 0), (4, 8), (6, 7)]),
            orbit(7, &[(5, 0)]),
        ],
        13 => vec![
            orbit(1, &[(2, 0)]),
            orbit(7, &[(3, 0)]),
            orbit(11, &[(4, 0)]),
            orbit(3, &[(3, 7), (3, 11), (3, 3)]),
        ],
        14 => vec![
            orbit(0, &[(3, 2)]),
            orbit(0, &[(8, 2)]),
            orbit(2, &[(4, 0), (1, 0)]),
            orbit(2, &[(5, 0), (3, 2)]),
            orbit(2, &[(6, 2), (4, 2)]),
            orbit(2, &[(8, 2), (7, 0), (5, 0)]),
            orbit(2, &[(6, 2)]),
            orbit(2, &[(2, 0), (6, 0)]),
        ],
        _ => return Err(PredefinedArQuiverError::Unsupported),
    };
    Ok(orbit_data(2, orbits))
}

fn bg1(r: usize) -> PredefinedArQuiverData {
    let mut data = rows![[4, 0], [3, 0]];
    for i in 1..=4 * r - 2 {
        data.push(vec![2 * i - 1, 2 * i + 4, 2 * i + 2]);
        data.push(vec![2 * i, 2 * i + 3, 2 * i + 1]);
    }
    data.push(vec![8 * r - 3, 8 * r]);
    data.push(vec![8 * r - 2, 8 * r - 1]);
    mesh_data(2, data)
}

fn bg3(r: usize) -> PredefinedArQuiverData {
    if r == 1 {
        return orbit_data(
            2,
            vec![
                orbit(1, &[(3, 0), (3, 2)]),
                orbit(1, &[(3, 0)]),
                orbit(3, &[(1, 1)]),
            ],
        );
    }
    let mut orbits = vec![
        orbit(r, &[(3, 0), (r + 2, 2)]),
        orbit(1, &[(r + 2, 0)]),
        orbit(r + 2, &[(4, r + 2)]),
    ];
    let mut orbit_number = 3;
    for _ in 4..=r + 1 {
        orbit_number += 1;
        orbits.push(orbit(
            r + 2,
            &[(orbit_number + 1, r + 2), (orbit_number - 1, 0)],
        ));
    }
    orbits.push(orbit(r + 2, &[(orbit_number, 0), (1, r)]));
    orbit_data(2, orbits)
}

fn bg4(r: usize) -> PredefinedArQuiverData {
    match r {
        1 => orbit_data(
            2,
            vec![
                orbit(4, &[(2, 0), (1, 2)]),
                orbit(2, &[(3, 0)]),
                orbit(3, &[(4, 1)]),
                orbit(1, &[(3, 0), (3, 2)]),
            ],
        ),
        2 => orbit_data(
            2,
            vec![
                orbit(5, &[(2, 0), (3, 0)]),
                orbit(3, &[(4, 0)]),
                orbit(2, &[(1, 2), (1, 5)]),
                orbit(4, &[(5, 4)]),
                orbit(4, &[(5, 2), (4, 0)]),
            ],
        ),
        3 => orbit_data(
            2,
            vec![
                orbit(6, &[(2, 0), (3, 3)]),
                orbit(4, &[(4, 0)]),
                orbit(6, &[(3, 3), (1, 3)]),
                orbit(5, &[(5, 5)]),
                orbit(5, &[(4, 0), (6, 2)]),
                orbit(2, &[(5, 0), (5, 3)]),
            ],
        ),
        _ if r % 2 == 0 => bg4_even(r),
        _ => bg4_odd(r),
    }
}

fn bg4_even(r: usize) -> PredefinedArQuiverData {
    let half_r_plus_two = usize::midpoint(r, 2);
    let mut orbits = vec![
        orbit(r + 3, &[(2, 0), (3, 3)]),
        orbit(r + 1, &[(3 + r / 2, 0)]),
        orbit(r + 3, &[(1, r), (4, 0)]),
    ];
    let mut orbit_number = 3;
    for _ in 4..=r / 2 + 1 {
        orbit_number += 1;
        orbits.push(orbit(
            r + 3,
            &[(orbit_number + 1, 0), (orbit_number - 1, r + 3)],
        ));
    }
    orbit_number += 1;
    orbits.push(orbit(
        half_r_plus_two,
        &[
            (orbit_number - 1, half_r_plus_two),
            (orbit_number - 1, r + 3),
        ],
    ));
    orbit_number += 1;
    orbits.push(orbit(r + 2, &[(orbit_number + 1, r + 2)]));
    for _ in 1..=half_r_plus_two - 2 {
        orbit_number += 1;
        orbits.push(orbit(
            r + 2,
            &[(orbit_number - 1, 0), (orbit_number + 1, r + 2)],
        ));
    }
    orbits.push(orbit(
        r + 2,
        &[(orbit_number, 0), (orbit_number + 1, half_r_plus_two)],
    ));
    orbit_data(2, orbits)
}

fn bg4_odd(r: usize) -> PredefinedArQuiverData {
    let half_r_plus_one = r.div_ceil(2);
    let half_r_plus_three = usize::midpoint(r, 3);
    let mut orbits = vec![
        orbit(r + 3, &[(2, 0), (3, 3)]),
        orbit(r + 1, &[((r - 1) / 2 + 3, 0)]),
        orbit(r + 3, &[(4, 0), (1, r)]),
    ];
    let mut orbit_number = 3;
    for _ in 1..=(r - 1) / 2 - 2 {
        orbit_number += 1;
        orbits.push(orbit(
            r + 3,
            &[(orbit_number - 1, r + 3), (orbit_number + 1, 0)],
        ));
    }
    orbit_number += 1;
    orbits.push(orbit(
        r + 3,
        &[(orbit_number - 1, r + 3), (orbit_number, half_r_plus_three)],
    ));
    orbit_number += 1;
    orbits.push(orbit(r + 2, &[(orbit_number + 1, r + 2)]));
    for _ in 1..=half_r_plus_three - 3 {
        orbit_number += 1;
        orbits.push(orbit(
            r + 2,
            &[(orbit_number - 1, 0), (orbit_number + 1, r + 2)],
        ));
    }
    orbit_number += 1;
    orbits.push(orbit(
        r + 2,
        &[(orbit_number - 1, 0), (orbit_number + 1, half_r_plus_one)],
    ));
    orbit_number += 1;
    orbits.push(orbit(
        half_r_plus_one,
        &[
            (orbit_number - 1, 0),
            (orbit_number - 1, half_r_plus_one + 1),
        ],
    ));
    orbit_data(2, orbits)
}

fn bg6(r: usize) -> PredefinedArQuiverData {
    if r == 1 {
        return orbit_data(
            2,
            vec![
                orbit(0, &[(4, 0)]),
                orbit(1, &[(3, 2)]),
                orbit(3, &[(4, 0), (2, 1)]),
                orbit(1, &[(3, 1), (3, 3)]),
            ],
        );
    }
    let mut orbits = vec![
        orbit(2 * r - 2, &[(4, 0)]),
        orbit(1, &[(3, 2 * r)]),
        orbit(2 * r + 1, &[(2, 1), (4, 2 * r)]),
    ];
    if r == 2 {
        orbits.push(orbit(2 * r + 1, &[(3, 1), (5, 2)]));
    } else {
        orbits.push(orbit(2 * r + 1, &[(3, 1), (5, 2 * r + 1)]));
    }
    let mut orbit_number = 4;
    for _ in 1..=r.saturating_sub(3) {
        orbit_number += 1;
        orbits.push(orbit(
            2 * r + 1,
            &[(orbit_number + 1, 2 * r + 1), (orbit_number - 1, 0)],
        ));
    }
    if r > 2 {
        orbit_number += 1;
        orbits.push(orbit(
            2 * r + 1,
            &[(orbit_number + 1, r), (orbit_number - 1, 0)],
        ));
    }
    orbits.push(orbit(r, &[(orbit_number, r + 1), (orbit_number, 0)]));
    orbit_data(2, orbits)
}

fn bg10(r: usize) -> PredefinedArQuiverData {
    let mut orbits = vec![
        orbit(0, &[(3, 0)]),
        orbit(0, &[(2 * r + 7, 1)]),
        orbit(1, &[(4, 1)]),
    ];
    let mut orbit_number = 3;
    for _ in 4..=2 * r + 6 {
        orbit_number += 1;
        orbits.push(orbit(1, &[(orbit_number + 1, 1), (orbit_number - 1, 0)]));
    }
    orbits.push(orbit(1, &[(2, 0), (2 * r + 6, 0), (2 * r + 8, 1)]));
    orbits.push(orbit(1, &[(2 * r + 7, 0)]));
    orbit_data(2, orbits)
}

fn bg2() -> PredefinedArQuiverData {
    mesh_data(
        2,
        rows![
            [11, 0],
            [19, 0],
            [12, 1],
            [13, 3],
            [14, 4],
            [15, 5],
            [16, 6],
            [17, 7],
            [18, 8],
            [20, 2],
            [21, 10],
            [22, 1, 11],
            [23, 3, 12],
            [24, 4, 13],
            [25, 5, 14],
            [26, 6, 15],
            [27, 7, 16],
            [28, 8, 17],
            [61, 66],
            [30, 2, 19],
            [31, 10, 20],
            [32, 11, 21],
            [33, 12, 22],
            [34, 13, 23],
            [35, 14, 24],
            [36, 15, 25],
            [37, 16, 26],
            [38, 17, 27],
            [39, 18, 28],
            [19, 54, 61],
            [40, 20, 30],
            [42, 21, 31],
            [44, 22, 32],
            [46, 23, 33],
            [48, 24, 34],
            [50, 25, 35],
            [52, 26, 36],
            [54, 27, 37],
            [40, 28, 38],
            [30, 55, 38, 54],
            [40, 55],
            [39, 31, 41, 40],
            [42, 41],
            [56, 32, 43, 42],
            [44, 43],
            [45, 57, 33, 44],
            [46, 45],
            [47, 58, 34, 46],
            [48, 47],
            [49, 59, 35, 48],
            [50, 49],
            [51, 60, 36, 50],
            [52, 51],
            [53, 61, 37, 52],
            [54, 53],
            [29, 42, 39],
            [62, 44, 56],
            [63, 46, 57],
            [64, 48, 58],
            [65, 50, 59],
            [66, 52, 60],
            [56, 29],
            [57, 62],
            [58, 63],
            [59, 64],
            [60, 65],
        ],
    )
}

fn bg5() -> PredefinedArQuiverData {
    mesh_data(
        2,
        rows![
            [11, 0],
            [17, 0],
            [12, 1],
            [13, 3],
            [14, 4],
            [18, 2],
            [19, 6],
            [20, 7],
            [21, 8],
            [22, 9],
            [23, 10],
            [24, 1, 11],
            [25, 3, 12],
            [26, 4, 13],
            [27, 5, 14],
            [63, 69],
            [28, 16],
            [29, 2, 17],
            [6, 30, 18],
            [7, 31, 19],
            [8, 32, 20],
            [9, 33, 21],
            [10, 34, 22],
            [11, 35, 23],
            [12, 36, 24],
            [13, 37, 25],
            [14, 38, 26],
            [47, 16, 63],
            [17, 39, 28],
            [18, 40, 29],
            [19, 41, 30],
            [20, 42, 31],
            [21, 43, 32],
            [22, 44, 33],
            [23, 45, 34],
            [24, 46, 35],
            [25, 47, 36],
            [26, 39, 37],
            [37, 48, 28, 47],
            [29, 49, 38, 39],
            [30, 50, 57, 40],
            [31, 51, 58, 41],
            [32, 52, 59, 42],
            [33, 53, 60, 43],
            [34, 54, 61, 44],
            [35, 55, 62, 45],
            [36, 56, 63, 46],
            [47, 56],
            [39, 48],
            [40, 49],
            [41, 50],
            [42, 51],
            [43, 52],
            [44, 53],
            [45, 54],
            [46, 55],
            [40, 27, 38],
            [41, 64, 57],
            [42, 65, 58],
            [43, 66, 59],
            [44, 67, 60],
            [45, 68, 61],
            [46, 69, 62],
            [57, 15, 27],
            [58, 70, 64],
            [59, 71, 65],
            [60, 72, 66],
            [61, 67],
            [62, 68],
            [64, 15],
            [65, 70],
            [66, 71],
        ],
    )
}

fn bg7() -> PredefinedArQuiverData {
    orbit_data(
        2,
        vec![
            orbit(12, &[(3, 0)]),
            orbit(1, &[(1, 7)]),
            orbit(15, &[(4, 0)]),
            orbit(15, &[(5, 0), (3, 15)]),
            orbit(15, &[(6, 0), (4, 15)]),
            orbit(15, &[(8, 0), (5, 15)]),
            orbit(15, &[(8, 15)]),
            orbit(15, &[(7, 0), (9, 0), (6, 15)]),
            orbit(15, &[(10, 0), (8, 15)]),
            orbit(15, &[(9, 15)]),
        ],
    )
}

fn bg8() -> PredefinedArQuiverData {
    let mut data = rows![
        [5, 0],
        [4, 16, 0],
        [6, 1],
        [7, 3],
        [2, 4],
        [1, 8, 5],
        [3, 9, 6],
        [5, 17, 2],
        [6, 18, 8],
        [7, 19, 9],
        [20, 10],
    ];
    for i in 1..=4 {
        data.push(vec![20 + i, 10 + i]);
    }
    data.extend(rows![[53, 57], [2, 26, 16], [8, 27, 17]]);
    for i in 1..=7 {
        data.push(vec![8 + i, 27 + i, 17 + i]);
    }
    data.extend(rows![[16, 41, 53], [17, 35, 26]]);
    for i in 1..=6 {
        data.push(vec![17 + i, 35 + i, 26 + i]);
    }
    data.extend(rows![
        [24, 35, 33],
        [26, 42, 33, 41],
        [34, 43, 27, 35],
        [49, 44, 28, 36],
    ]);
    for i in 1..=4 {
        data.push(vec![49 + i, 44 + i, 28 + i, 36 + i]);
    }
    data.extend(rows![[41, 48], [35, 42]]);
    for i in 1..=5 {
        data.push(vec![35 + i, 42 + i]);
    }
    data.extend(rows![[25, 36, 34], [54, 37, 49]]);
    for i in 1..=3 {
        data.push(vec![54 + i, 37 + i, 49 + i]);
    }
    data.extend(rows![[49, 25], [50, 54], [51, 55], [52, 56]]);
    mesh_data(2, data)
}

fn bg9() -> PredefinedArQuiverData {
    mesh_data(
        2,
        rows![
            [5, 0],
            [21, 26, 0],
            [6, 1],
            [11, 7],
            [8, 4],
            [1, 9, 5],
            [3, 10, 6],
            [4, 12, 19, 11],
            [5, 13, 16, 8],
            [6, 14, 17, 9],
            [7, 15, 18, 10],
            [11, 15],
            [8, 12],
            [9, 13],
            [10, 14],
            [8, 23, 19],
            [9, 20, 16],
            [10, 21, 17],
            [11, 22, 18],
            [16, 25, 23],
            [17, 20],
            [2, 18, 21],
            [19, 24, 22],
            [22, 27, 2],
            [23, 28, 24],
            [25, 28],
            [2, 26],
            [24, 27],
        ],
    )
}

fn bg11() -> PredefinedArQuiverData {
    let mut data = rows![[14, 0], [10, 0], [11, 2], [15, 1], [16, 4]];
    for i in 1..=5 {
        data.push(vec![16 + i, 4 + i]);
    }
    data.extend(rows![
        [2, 22, 10],
        [3, 23, 11],
        [55, 60],
        [24, 13],
        [1, 25, 14],
        [4, 26, 15],
    ]);
    for i in 1..=7 {
        data.push(vec![4 + i, 26 + i, 15 + i]);
    }
    data.extend(rows![[13, 41, 55], [14, 34, 24]]);
    for i in 1..=7 {
        data.push(vec![14 + i, 34 + i, 24 + i]);
    }
    data.extend(rows![
        [22, 34, 32],
        [24, 32, 42, 41],
        [25, 33, 43, 34],
        [26, 44, 50, 35],
    ]);
    for i in 1..=5 {
        data.push(vec![26 + i, 44 + i, 50 + i, 35 + i]);
    }
    data.extend(rows![[41, 49], [34, 42]]);
    for i in 1..=6 {
        data.push(vec![34 + i, 42 + i]);
    }
    data.extend(rows![[23, 35, 33], [56, 36, 50]]);
    for i in 1..=4 {
        data.push(vec![56 + i, 36 + i, 50 + i]);
    }
    data.extend(rows![
        [12, 50, 23],
        [51, 61, 56],
        [52, 62, 57],
        [53, 63, 58],
        [54, 59],
        [56, 12],
        [57, 61],
        [58, 62],
    ]);
    mesh_data(2, data)
}

fn bg12() -> PredefinedArQuiverData {
    let mut data = rows![
        [21, 0],
        [6, 0],
        [7, 2],
        [8, 3],
        [22, 1],
        [23, 5],
        [2, 24, 6],
        [3, 25, 7],
        [4, 26, 8],
        [27, 9],
    ];
    for i in 1..=7 {
        data.push(vec![27 + i, 9 + i]);
    }
    data.extend(rows![
        [17, 35, 34],
        [36, 18],
        [37, 19],
        [38, 20],
        [1, 39, 21],
        [5, 40, 22],
    ]);
    for i in 1..=11 {
        data.push(vec![5 + i, 40 + i, 22 + i]);
    }
    data.extend(rows![[34, 52, 51], [18, 53, 35]]);
    for i in 1..=15 {
        data.push(vec![18 + i, 53 + i, 35 + i]);
    }
    data.extend(rows![[51, 69, 68], [35, 70, 52]]);
    for i in 1..=15 {
        data.push(vec![35 + i, 70 + i, 52 + i]);
    }
    data.extend(rows![[68, 103, 85], [52, 104, 69]]);
    for i in 1..=15 {
        data.push(vec![52 + i, 104 + i, 69 + i]);
    }
    data.extend(rows![[119, 102], [103, 86]]);
    for i in 1..=15 {
        data.push(vec![103 + i, 86 + i]);
    }
    data.extend(rows![[85, 86, 120, 119], [69, 87, 121, 103]]);
    for i in 1..=15 {
        data.push(vec![69 + i, 87 + i, 121 + i, 103 + i]);
    }
    data.extend(rows![[119, 137, 136], [103, 138, 120]]);
    for i in 1..=15 {
        data.push(vec![103 + i, 138 + i, 120 + i]);
    }
    data.extend(rows![[136, 153], [120, 137]]);
    for i in 1..=15 {
        data.push(vec![120 + i, 137 + i]);
    }
    mesh_data(2, data)
}

fn bg13() -> PredefinedArQuiverData {
    mesh_data(
        2,
        rows![
            [2, 0],
            [8, 0],
            [4, 1],
            [1, 9, 2],
            [3, 10, 4],
            [11, 5],
            [12, 6],
            [13, 20],
            [2, 15, 8],
            [4, 17, 9],
            [5, 19, 10],
            [6, 13, 11],
            [11, 20, 23, 19],
            [13, 26, 23],
            [8, 12, 14, 13],
            [15, 14],
            [9, 16, 21, 15],
            [17, 16],
            [10, 18, 22, 17],
            [19, 18],
            [7, 15, 12],
            [17, 24, 21],
            [19, 25, 22],
            [21, 7],
            [22, 24],
            [23, 25],
        ],
    )
}

fn bg14() -> PredefinedArQuiverData {
    mesh_data(
        2,
        rows![
            [5, 0],
            [20, 0],
            [1, 6, 5],
            [7, 3],
            [8, 4],
            [5, 9, 8],
            [3, 10, 6],
            [4, 11, 7],
            [8, 14, 11],
            [6, 12, 9],
            [7, 13, 10],
            [9, 15, 20, 14],
            [10, 16, 18, 12],
            [11, 17, 19, 13],
            [14, 17],
            [12, 15],
            [13, 16],
            [2, 12, 20],
            [13, 18],
            [14, 19],
        ],
    )
}

fn d4_subspace() -> PredefinedArQuiverData {
    mesh_data(
        4,
        rows![
            [0],
            [1, 0],
            [1, 0],
            [1, 0],
            [2, 3, 4, 1],
            [5, 2],
            [5, 3],
            [5, 4],
            [6, 7, 8, 5],
            [9, 6],
            [9, 7],
            [9, 8],
        ],
    )
}

fn e6_subspace() -> PredefinedArQuiverData {
    mesh_data(
        6,
        rows![
            [0],
            [1, 0],
            [2, 0],
            [1, 0],
            [1, 0],
            [5, 0],
            [2, 4, 5, 1],
            [3, 7, 2],
            [8, 3],
            [7, 4],
            [7, 6, 5],
            [11, 6],
            [8, 10, 11, 7],
            [9, 13, 8],
            [14, 9],
            [13, 10],
            [13, 12, 11],
            [17, 12],
            [14, 16, 17, 13],
            [15, 19, 14],
            [20, 15],
            [19, 16],
            [19, 18, 17],
            [23, 18],
            [20, 22, 23, 19],
            [21, 25, 20],
            [26, 21],
            [25, 22],
            [25, 24, 23],
            [29, 24],
            [26, 28, 29, 25],
            [27, 31, 26],
            [32, 27],
            [31, 28],
            [31, 30, 29],
            [35, 30],
        ],
    )
}

fn a3_zero() -> PredefinedArQuiverData {
    mesh_data(3, rows![[0], [1, 0], [4, 0], [2, 1], [3, 4]])
}

fn a12_zero() -> PredefinedArQuiverData {
    mesh_data(
        3,
        rows![
            [0],
            [1, 0],
            [4, 1, 0],
            [7, 5],
            [3, 4],
            [3, 2, 1],
            [5, 6, 3],
            [6, 2],
            [7, 8, 6],
        ],
    )
}

fn a22_commutative() -> PredefinedArQuiverData {
    mesh_data(
        4,
        rows![
            [0],
            [1, 0],
            [1, 0],
            [5, 0],
            [2, 3, 1],
            [5, 2],
            [5, 3],
            [6, 7, 4, 5],
            [8, 6],
            [8, 7],
            [9, 10, 8],
        ],
    )
}

fn r_nilpotent() -> PredefinedArQuiverData {
    mesh_data(
        2,
        rows![
            [4, 0],
            [1, 0],
            [1, 6, 4],
            [3, 6],
            [3, 2, 1],
            [4, 5, 3],
            [5, 2],
        ],
    )
}

#[cfg(test)]
mod tests {
    use serde_json::{Value, json};

    use super::{
        PREDEFINED_AR_QUIVER_CATALOG, PredefinedArQuiver, PredefinedArQuiverData,
        PredefinedArQuiverError, predefined_ar_quiver, predefined_ar_quiver_data,
    };
    use crate::ArQuiverDescription;

    fn cases() -> Vec<(&'static str, PredefinedArQuiver)> {
        use PredefinedArQuiver::{
            A3Zero, A12Zero, A22Commutative, BongartzGabriel, BongartzGabrielOrbits, D4Subspace,
            E6Subspace, RNilpotent,
        };

        let mut cases = vec![
            ("D4 subspace", D4Subspace),
            ("E6 subspace", E6Subspace),
            ("A3 zero", A3Zero),
            ("A1,2 zero", A12Zero),
            ("A2,2 comm", A22Commutative),
            ("R nilp", RNilpotent),
        ];
        for number in [2, 5, 7, 8, 9, 11, 12, 13, 14] {
            cases.push((
                match number {
                    2 => "BG2",
                    5 => "BG5",
                    7 => "BG7",
                    8 => "BG8",
                    9 => "BG9",
                    11 => "BG11",
                    12 => "BG12",
                    13 => "BG13",
                    14 => "BG14",
                    _ => unreachable!(),
                },
                BongartzGabriel {
                    number,
                    parameter: None,
                },
            ));
        }
        for number in [2, 13, 14] {
            cases.push((
                match number {
                    2 => "BGo2",
                    13 => "BGo13",
                    14 => "BGo14",
                    _ => unreachable!(),
                },
                BongartzGabrielOrbits { number },
            ));
        }
        for (label, number, parameter) in [
            ("BG1r1", 1, 1),
            ("BG1r2", 1, 2),
            ("BG3r1", 3, 1),
            ("BG3r2", 3, 2),
            ("BG3r3", 3, 3),
            ("BG4r1", 4, 1),
            ("BG4r2", 4, 2),
            ("BG4r3", 4, 3),
            ("BG4r4", 4, 4),
            ("BG4r5", 4, 5),
            ("BG4r6", 4, 6),
            ("BG4r7", 4, 7),
            ("BG6r1", 6, 1),
            ("BG6r2", 6, 2),
            ("BG6r3", 6, 3),
            ("BG6r4", 6, 4),
            ("BG10r1", 10, 1),
            ("BG10r2", 10, 2),
        ] {
            cases.push((
                label,
                BongartzGabriel {
                    number,
                    parameter: Some(parameter),
                },
            ));
        }
        cases
    }

    fn qpa_data(data: PredefinedArQuiverData) -> Value {
        let description = match data.description {
            ArQuiverDescription::Meshes(meshes) => Value::Array(
                meshes
                    .into_iter()
                    .map(|mesh| {
                        let mut row = mesh
                            .direct_predecessors
                            .into_iter()
                            .map(|position| json!(position + 1))
                            .collect::<Vec<_>>();
                        row.push(json!(mesh.tau.map_or(0, |position| position + 1)));
                        Value::Array(row)
                    })
                    .collect(),
            ),
            ArQuiverDescription::Orbits(orbits) => json!([
                "orbits",
                orbits
                    .into_iter()
                    .map(|orbit| {
                        let mut row = vec![json!(orbit.length_minus_one)];
                        row.extend(
                            orbit
                                .direct_predecessors
                                .into_iter()
                                .map(|position| json!([position.orbit + 1, position.shift])),
                        );
                        row
                    })
                    .collect::<Vec<_>>()
            ]),
        };
        json!([
            data.number_of_indecomposables,
            data.number_of_projectives,
            description
        ])
    }

    #[test]
    fn complete_pinned_qpa_predefined_catalog_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/degorderfinitetype_predefined.schema-v1.json"
        ))
        .unwrap();
        let actual = cases()
            .into_iter()
            .map(|(label, selector)| {
                json!([
                    label,
                    qpa_data(predefined_ar_quiver_data(selector).unwrap())
                ])
            })
            .collect::<Vec<_>>();

        assert_eq!(json!(actual), expected["cases"]);
    }

    #[test]
    fn catalog_text_matches_pinned_qpa_output() {
        assert_eq!(
            PREDEFINED_AR_QUIVER_CATALOG,
            include_str!("../../../tests/fixtures/degorderfinitetype_predefined_catalog.txt")
        );
    }

    #[test]
    fn every_oracle_catalog_case_constructs_a_numerical_ar_quiver() {
        for (label, selector) in cases() {
            predefined_ar_quiver(selector)
                .unwrap_or_else(|error| panic!("QPA catalog case {label} failed: {error}"));
        }
    }

    #[test]
    fn qpa_dispatch_domain_and_constructor_are_preserved() {
        use PredefinedArQuiver::{BongartzGabriel, BongartzGabrielOrbits};

        for selector in [
            BongartzGabriel {
                number: 1,
                parameter: None,
            },
            BongartzGabriel {
                number: 2,
                parameter: Some(1),
            },
            BongartzGabriel {
                number: 4,
                parameter: Some(0),
            },
            BongartzGabrielOrbits { number: 7 },
        ] {
            assert_eq!(
                predefined_ar_quiver_data(selector),
                Err(PredefinedArQuiverError::Unsupported)
            );
        }

        let ar = predefined_ar_quiver(PredefinedArQuiver::D4Subspace).unwrap();
        assert_eq!(ar.number_of_indecomposables(), 12);
        assert_eq!(ar.number_of_projectives(), 4);
    }
}
