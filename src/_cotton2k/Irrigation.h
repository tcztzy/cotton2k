#ifndef IRRIGATION_TYPE
#define IRRIGATION_TYPE
typedef struct Irrigation
{
    int day;                // date of application (DOY)
    int method;             // index of irrigation method ( 0 = sprinkler; 1 = furrow; 2 = drip)
    int LocationColumnDrip; // horizontal placement of side-dressed fertilizer, cm.
    int LocationLayerDrip;  // vertical placement of side-dressed fertilizer, cm.
    double amount;          // water applied, mm.
} Irrigation;
#endif
