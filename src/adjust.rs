use crate::{
    nadj, AdjAddHeightRate, AdjAddMSNodesRate, AdjAddSitesRate, AdjGreenBollAbsc, AdjSquareAbsc,
    DayOfSimulation, Daynum, FirstSquare, GoBack, MapDataAllSiteNum, MapDataGreenBollNum,
    MapDataMainStemNodes, MapDataPlantHeight, MapDataSquareNum, NumAdjustDays, NumFruitSites,
    NumGreenBolls, NumSquares, PlantHeight, Scratch21,
};

/// This function adjusts plant height and plant fruiting map, when data for such adjustments are available.
///
/// This function is called from [Profile::adjust]. it calls [GoBack()].
///
/// The following global variables are referenced here:
/// * [Daynum]
/// * [MapDataAllSiteNum]
/// * [MapDataGreenBollNum]
/// * [MapDataMainStemNodes]
/// * [MapDataPlantHeight]
/// * [MapDataSquareNum]
/// * [NumAdjustDays]
/// * [NumFruitSites]
/// * [NumGreenBolls]
/// * [NumSquares]
/// * [PlantHeight]
///
/// The following global variables are set:
//  AdjAddHeightRate, AdjAddMSNodesRate, AdjAddSitesRate, AdjGreenBollAbsc,
//  AdjSquareAbsc, FirstSquare, nadj.
//     The following arguments are used:
//  i - the number of this adjustment, as read from the *.MAP file.
//  jj - the type of this adjustment.
pub unsafe fn PlantAdjustments(i: usize, jj: i32) {
    // Define nadj(jj) as true or false, where jj is:
    //  0 - for main stem nodes.
    //  1 - for height.
    //  2 - for total site number.
    //  3 - for square number.
    //  4 - for green boll number.
    //  5 - for open boll number.
    match jj {
        0 => {
            // adjust the number of main stem nodes
            if MapDataMainStemNodes[i] <= 0. {
                return;
            }
            // Compute targeted number of adjusted number of fruiting branches (mntarget).
            // The difference between the measured and the simulated number of nodes should be more than one node. Otherwise no adjustments are made.

            // simulated number of main stem node count.
            let mnsim = Scratch21[(DayOfSimulation - 1) as usize].mainStemNodes as f64;
            if (MapDataMainStemNodes[i] - mnsim).abs() <= 1. {
                nadj[0] = false;
                return;
            }
            nadj[0] = true;
            let mntarget = MapDataMainStemNodes[i]; //  targeted number of main stem nodes.

            // Compute the ratio used to adjust the rate of formation of mainstem nodes.

            // simulated number of main stem nodes, at the start of plant adjustment period.
            let mn0 =
                Scratch21[(DayOfSimulation - NumAdjustDays - 1) as usize].mainStemNodes as f64;
            AdjAddMSNodesRate = 1.;
            if mnsim != mn0 {
                AdjAddMSNodesRate = (mntarget - mn0) / (mnsim - mn0) as f64;
                if AdjAddMSNodesRate > 0.98 && AdjAddMSNodesRate < 1.02 {
                    AdjAddMSNodesRate = 1.;
                }
                if AdjAddMSNodesRate < 0. {
                    AdjAddMSNodesRate = 0.;
                }
            }
            if AdjAddMSNodesRate == 1. {
                nadj[0] = false;
            } else {
                nadj[0] = true;
                // AdjAddMSNodesRate will be used in function AddFruitingBranch()
                println!(
                    " Apply plant adjustment for main stem nodes to date {}",
                    Daynum,
                );
                GoBack();
            }
        }
        1 => {
            // Plant stem height
            if MapDataPlantHeight[i] <= 0. {
                return;
            }
            // Compute targeted adjusted plant height, if difference is larger than 5% of the average.
            // Note that adjustment is by 90% of the difference.
            // simulated plant height.
            let zsim = PlantHeight;
            // the difference between simulated and plant adjustment values.
            let ddif = MapDataPlantHeight[i] - PlantHeight;
            // The difference should be at least 5% of the height, otherwise no adjustments made.
            if ddif.abs() <= (MapDataPlantHeight[i] + PlantHeight) / 40. {
                nadj[1] = false;
                return;
            }
            // targeted plant height.
            let ztarget = PlantHeight + 0.9 * ddif;
            nadj[1] = true;
            // Compute the ratio to adjust rate of growth of plant height.
            AdjAddHeightRate = 1.;
            // plant height at start of adjustment
            let pHeight = Scratch21[(DayOfSimulation - NumAdjustDays - 1) as usize].plantHeight;
            if (zsim - pHeight).abs() > 0. {
                AdjAddHeightRate = (ztarget - pHeight) / (zsim - pHeight);
            }
            // no negative growth rates for plant height.
            if AdjAddHeightRate < 0. {
                AdjAddHeightRate = 0.;
            }
            // If AdjAddHeightRate value is near 1, no adjustments are made.
            if AdjAddHeightRate > 0.98 && AdjAddHeightRate < 1.02 {
                nadj[1] = false;
            } else {
                nadj[1] = true;
                // AdjAddHeightRate will be used in function AddPlantHeight()
                println!(" Apply plant adjustment for stem height to date {}", Daynum,);
                GoBack();
            }
        }
        2 => {
            // total number of fruiting sites
            if MapDataAllSiteNum[i] <= 0. {
                return;
            }
            // If first square has not yet been simulated, but map data shows squares,
            // adjust the day of first square.
            if FirstSquare <= 0 {
                FirstSquare = (Daynum as f64 - 3. * MapDataAllSiteNum[i]) as i32;
                nadj[2] = true;
                return;
            }
            // simulated number of total fruiting sites.
            let sitesim = NumFruitSites;
            // Compute targeted adjusted total site number.
            if (MapDataAllSiteNum[i] - NumFruitSites as f64).abs() <= 2. {
                nadj[2] = false;
                return;
            }
            nadj[2] = true;
            // targeted number of total fruiting sites.
            let sitarget = MapDataAllSiteNum[i];
            // Compute the ratio to adjust the rate of formation of sites on fruiting branches.
            if sitarget <= 0. {
                AdjAddSitesRate = 1.;
            }
            // number of sites at start of adjustment
            let pNumFruitSites =
                Scratch21[(DayOfSimulation - NumAdjustDays - 1) as usize].numFruitSites;
            if sitesim != pNumFruitSites {
                AdjAddSitesRate =
                    (sitarget - pNumFruitSites as f64) / (sitesim - pNumFruitSites) as f64;
                if AdjAddSitesRate > 0.98 && AdjAddSitesRate < 1.02 {
                    AdjAddSitesRate = 1.;
                }
                if AdjAddSitesRate < 0. {
                    AdjAddSitesRate = 0.;
                }
            }
            if AdjAddSitesRate == 1. {
                nadj[2] = false;
            } else {
                nadj[2] = true;
                // AdjAddSitesRate will be used in function AddFruitingNode()
                println!(
                    " Apply plant adjustment for total number of sites to date {}",
                    Daynum
                );
                GoBack();
            }
        }
        3 => {
            // number of squares
            if MapDataSquareNum[i] <= 0. {
                return;
            }
            // Adjust actual square numbers. Amount to add each day is the difference between target and actual simulated value on adjustment date divided by the length of the  adjustment period.
            if (NumSquares - MapDataSquareNum[i]) > 1. && NumAdjustDays > 0 {
                nadj[3] = true;
                AdjSquareAbsc =
                    1. - (MapDataSquareNum[i] / NumSquares).powf(1. / NumAdjustDays as f64);
                //     AdjSquareAbsc will be used in function AdjustAbscission()
                println!(
                    " Apply plant adjustment for number of squares to date {}",
                    Daynum
                );
                GoBack();
            } else {
                nadj[3] = false;
            }
        }
        4 => {
            // number of green bolls
            if MapDataGreenBollNum[i] <= 0. {
                return;
            }
            // Adjust actual green boll numbers using the same method as for squares.
            if (NumGreenBolls - MapDataGreenBollNum[i]) > 1. && NumAdjustDays > 0 {
                AdjGreenBollAbsc =
                    1. - (MapDataGreenBollNum[i] / NumGreenBolls).powf(1. / NumAdjustDays as f64);
                nadj[4] = true;
                // AdjGreenBollAbsc will be used in function AdjustAbscission()
                println!(
                    " Apply plant adjustment for number of green bolls to date {}",
                    Daynum,
                );
                GoBack();
            } else {
                nadj[4] = false;
            }
        }
        _ => {}
    }
}
