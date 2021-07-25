from math import exp


# pylint: disable=no-member, too-few-public-methods
class Photosynthesis:
    # pylint: disable=too-many-arguments
    def column_shading(
        self,
        row_space,
        plant_row_column,
        column_width,
        max_leaf_area_index,
        relative_radiation_received_by_a_soil_column,
    ):
        zint = 1.0756 * self.plant_height / row_space
        for k in range(20):
            sw = (k + 1) * column_width
            if k <= plant_row_column:
                k0 = plant_row_column - k
                sw1 = sw - column_width / 2
            else:
                sw1 = sw - column_width / 2 - (plant_row_column + 1) * column_width
                k0 = k
            shade = 0
            if sw1 < self.plant_height:
                shade = 1 - (sw1 / self.plant_height) ** 2
                if (
                    self.light_interception < zint
                    and self.leaf_area_index < max_leaf_area_index
                ):
                    shade *= self.light_interception / zint
            relative_radiation_received_by_a_soil_column[k0] = max(0.05, 1 - shade)

    def compute_light_interception(
        self,
        max_leaf_area_index: float,
        row_space: float,
    ):
        if self.version < 0x0500:  # type: ignore[attr-defined]
            zint = 1.0756 * self.plant_height / row_space  # type: ignore[attr-defined]
            lfint = (
                0.80 * self.leaf_area_index  # type: ignore[attr-defined]
                if self.leaf_area_index <= 0.5  # type: ignore[attr-defined]
                else 1 - exp(0.07 - 1.16 * self.leaf_area_index)  # type: ignore
            )
            if lfint > zint:
                light_interception = (zint + lfint) / 2
            elif self.leaf_area_index < max_leaf_area_index:  # type: ignore
                light_interception = lfint
            else:
                light_interception = zint
            return light_interception if light_interception < 1 else 1
        param = max(1.16, -0.1 * self.plant_height + 8)  # type: ignore[attr-defined]
        return 1 - exp(-param * self.leaf_area_index)  # type: ignore[attr-defined]
