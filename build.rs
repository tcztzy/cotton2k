extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    let cpp_sources = vec![
        "CottonPhenology.cpp",
        "FruitAbscission.cpp",
        "GeneralFunctions.cpp",
        "GettingInput_2.cpp",
        "global.cpp",
        "LeafAbscission.cpp",
        "PlantAdjustment.cpp",
        "PlantGrowth_2.cpp",
        "PlantGrowth_3.cpp",
        "PlantNitrogen.cpp",
        "RootGrowth_1.cpp",
        "RootGrowth_2.cpp",
        "SoilNitrogen.cpp",
        "SoilProcedures_1.cpp",
        "SoilProcedures_2.cpp",
        "SoilProcedures_3.cpp",
        "SoilTemperature_2.cpp",
        "SoilTemperature_3.cpp",
    ];

    let x = cpp_sources.clone();
    cc::Build::new()
        .cpp(true)
        .files(cpp_sources)
        .compile("cotton2k");
    println!("cargo:rustc-link-lib=cotton2k");

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    for &s in x.iter() {
        println!("cargo:rerun-if-changed={}", s);
    }
    println!("cargo:rerun-if-changed=global.h");
    println!("cargo:rerun-if-changed=CottonSimulation.h");
    println!("cargo:rerun-if-changed=GeneralFunctions.h");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .clang_arg("-xc++")
        .clang_arg("-std=c++14")
        .header("global.h")
        .header("CottonSimulation.h")
        .header("GeneralFunctions.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
