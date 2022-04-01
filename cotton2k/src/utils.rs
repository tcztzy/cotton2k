use crate::Cotton2KError;

fn depth_of_layer(l: usize) -> Result<f64, Cotton2KError> {
    if l >= 40 {
        Err(Cotton2KError {
            level: 1,
            message: String::from("Out of index of soil layers!"),
        })
    } else {
        Ok((match l {
            0 => 2,
            1 => 2,
            2 => 2,
            3 => 4,
            38 => 10,
            39 => 10,
            _ => 5,
        }) as f64)
    }
}

fn width_of_column(k: usize, row_space: f64) -> Result<f64, Cotton2KError> {
    if k >= 20 {
        Err(Cotton2KError {
            level: 1,
            message: String::from("Out of index of soil columns!"),
        })
    } else {
        Ok(row_space / 20.)
    }
}

fn depth(l: usize) -> Result<f64, Cotton2KError> {
    if l >= 40 {
        Err(Cotton2KError {
            level: 1,
            message: String::from("Out of index of soil layers!"),
        })
    } else {
        Ok(match l {
            0 => 1.,
            1 => 3.,
            2 => 5.,
            3 => 8.,
            38 => 185.,
            39 => 195.,
            _ => (l * 5) as f64 - 7.5,
        })
    }
}

#[test]
fn test_depth() -> Result<(), Cotton2KError> {
    assert_eq!(depth(0)?, 1.);
    assert_eq!(depth(1)?, 3.);
    assert_eq!(depth(2)?, 5.);
    assert_eq!(depth(3)?, 8.);
    assert_eq!(depth(4)?, 12.5);
    assert_eq!(depth(5)?, 17.5);
    assert_eq!(depth(36)?, 172.5);
    assert_eq!(depth(37)?, 177.5);
    assert_eq!(depth(38)?, 185.);
    assert_eq!(depth(39)?, 195.);
    Ok(())
}

fn width(k: usize, row_space: f64) -> Result<f64, Cotton2KError> {
    if k >= 20 {
        Err(Cotton2KError {
            level: 1,
            message: String::from("Out of index of soil columns!"),
        })
    } else if row_space < 0. {
        Err(Cotton2KError {
            level: 1,
            message: String::from("Negative row space!"),
        })
    } else {
        Ok((k as f64 + 0.5) * row_space / 20.)
    }
}

///  It is called from [crate::DripFlow()].
pub fn cell_distance(
    l: usize,
    k: usize,
    l0: usize,
    k0: usize,
    row_space: f64,
) -> Result<f64, Cotton2KError> {
    Ok(
        ((depth(l)? - depth(l0)?).powi(2) + (width(k, row_space)? - width(k0, row_space)?).powi(2))
            .sqrt(),
    )
}

#[test]
fn test_cell_distance() -> Result<(), Cotton2KError> {
    assert_eq!(cell_distance(0, 0, 0, 0, 100.)?, 0f64);
    assert_eq!(cell_distance(0, 0, 0, 19, 100.)?, 95f64);
    assert_eq!(cell_distance(0, 0, 39, 0, 100.)?, 194f64);
    assert_eq!(
        cell_distance(0, 0, 39, 19, 100.)?,
        (194f64.powi(2) + 95f64.powi(2)).sqrt()
    );
    Ok(())
}

pub fn slab_vertical_location(distance: f64) -> Result<usize, Cotton2KError> {
    if distance > 200. || distance < 0. {
        Err(Cotton2KError {
            level: 1,
            message: String::from("distance out of boundary!"),
        })
    } else {
        let mut result: usize = 0;

        for l in 0..40 {
            if (depth(l)? + depth_of_layer(l)? / 2.) >= distance {
                result = l;
                break;
            }
        }

        Ok(result)
    }
}

#[test]
fn test_slab_vertical_location() -> Result<(), Cotton2KError> {
    assert_eq!(slab_vertical_location(5.)?, 2);
    assert_eq!(slab_vertical_location(10.)?, 3);
    assert_eq!(slab_vertical_location(108.)?, 23);
    assert_eq!(slab_vertical_location(200.)?, 39);
    Ok(())
}

pub fn slab_horizontal_location(distance: f64, row_space: f64) -> Result<usize, Cotton2KError> {
    if distance > row_space || distance < 0. {
        Err(Cotton2KError {
            level: 1,
            message: String::from("distance out of boundary!"),
        })
    } else {
        let mut result: usize = 0;

        for l in 0..20 {
            if (width(l, row_space)? + width_of_column(l, row_space)? / 2.) >= distance {
                result = l;
                break;
            }
        }

        Ok(result)
    }
}
#[test]
fn test_slab_horizontal_location() -> Result<(), Cotton2KError> {
    assert_eq!(slab_horizontal_location(10., 100.)?, 1);
    assert_eq!(slab_horizontal_location(20., 100.)?, 3);
    Ok(())
}

/// parameters used to correct photosynthesis for ambient CO2 concentration.
const CO2_PARAMETER: [f64; 45] = [
    1.0235, 1.0264, 1.0285, 1.0321, 1.0335, 1.0353, 1.0385, 1.0403, 1.0431, 1.0485, 1.0538, 1.0595,
    1.0627, 1.0663, 1.0716, 1.0752, 1.0784, 1.0823, 1.0880, 1.0923, 1.0968, 1.1019, 1.1087, 1.1172,
    1.1208, 1.1243, 1.1311, 1.1379, 1.1435, 1.1490, 1.1545, 1.1601, 1.1656, 1.1712, 1.1767, 1.1823,
    1.1878, 1.1934, 1.1990, 1.2045, 1.2101, 1.2156, 1.2212, 1.2267, 1.2323,
];
pub fn ambient_CO2_factor(year: i32) -> f64 {
    let co2indx = year - 1960; // count of years from 1960.
    if co2indx < 0 {
        1.
    } else if (co2indx as usize) < CO2_PARAMETER.len() {
        // for years 1960 to 2004
        CO2_PARAMETER[co2indx as usize]
    } else {
        // extrapolate for years after 2004
        CO2_PARAMETER.last().unwrap() + 0.004864 * (co2indx - 45) as f64
    }
}

#[test]
fn test_ambient_CO2_factor() {
    assert_eq!(ambient_CO2_factor(1959), 1.);
    assert_eq!(ambient_CO2_factor(1960), CO2_PARAMETER[0]);
}
