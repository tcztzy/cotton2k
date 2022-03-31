use crate::Cotton2KError;

fn depth_of_layer(l: usize) -> Result<f64, Cotton2KError> {
    if l >= 40 {
        Err(Cotton2KError {
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
            message: String::from("Out of index of soil columns!"),
        })
    } else {
        Ok(row_space / 20.)
    }
}

fn depth(l: usize) -> Result<f64, Cotton2KError> {
    if l >= 40 {
        Err(Cotton2KError {
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
            message: String::from("Out of index of soil columns!"),
        })
    } else if row_space < 0. {
        Err(Cotton2KError {
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
