use super::{Stage_GreenBoll, Stage_MatureBoll, Stage_Square, Stage_YoungGreenBoll, State};

#[no_mangle]
extern "C" fn ComputeSiteNumbers(state: &mut State)
//     This function calculates square, green boll, open boll, and abscised site numbers
//  (NumSquares, NumGreenBolls, NumOpenBolls, and AbscisedFruitSites, respectively), as
//  the sums of FruitFraction in all sites with appropriate FruitingCode.
//  It is called from function FruitingSitesAbscission().
{
    state.number_of_squares = 0f64;
    state.number_of_green_bolls = 0f64;
    state.number_of_open_bolls = 0f64;
    for k in 0..state.number_of_vegetative_branches as usize {
        for l in 0..state.vegetative_branches[k].number_of_fruiting_branches as usize {
            for m in 0..state.vegetative_branches[k].fruiting_branches[l].number_of_fruiting_nodes
                as usize
            {
                let site = &mut state.vegetative_branches[k].fruiting_branches[l].nodes[m];
                if site.stage == Stage_Square {
                    state.number_of_squares += site.fraction;
                } else if site.stage == Stage_YoungGreenBoll || site.stage == Stage_GreenBoll {
                    state.number_of_green_bolls += site.fraction;
                } else if site.stage == Stage_MatureBoll {
                    state.number_of_open_bolls += site.fraction;
                }
            }
        }
    }
    state.abscised_fruit_sites = state.number_of_fruiting_sites as f64
        - state.number_of_squares
        - state.number_of_green_bolls
        - state.number_of_open_bolls;
}
