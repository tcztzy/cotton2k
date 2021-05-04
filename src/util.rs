use chrono::prelude::*;
use chrono::Duration;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

/// This function sorts an array of values by its value (larger is first)
/// together with three indexes associated with each value.
#[no_mangle]
extern "C" fn SortArray(size: usize, data: *mut f64, ik: *mut i32, il: *mut i32, im: *mut i32) {
    let _data = unsafe { slice::from_raw_parts_mut(data, size) };
    let _ik = unsafe { slice::from_raw_parts_mut(ik, size) };
    let _il = unsafe { slice::from_raw_parts_mut(il, size) };
    let _im = unsafe { slice::from_raw_parts_mut(im, size) };
    let mut x = [(0f64, 0i32, 0i32, 0i32)].repeat(size);
    for i in 0..size {
        x[i] = (_data[i], _ik[i], _il[i], _im[i]);
    }
    x.sort_by(|a, b| (b.0).partial_cmp(&a.0).unwrap());
    unsafe {
        for (i, _x) in x.iter().enumerate() {
            *(data.add(i)) = _x.0;
            *(ik.add(i)) = _x.1;
            *(il.add(i)) = _x.2;
            *(im.add(i)) = _x.3;
        }
    }
}
