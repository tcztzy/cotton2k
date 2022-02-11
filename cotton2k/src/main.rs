use cotton2k::run;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        panic!("profile file path should be provided!");
    }
    run(std::path::Path::new(&args[1]));
}