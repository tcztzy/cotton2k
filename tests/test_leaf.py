from _cotton2k.leaf import temperature_on_leaf_growth_rate


def test_temperature_on_leaf_growth_rate():
    assert temperature_on_leaf_growth_rate(12) == 0
    assert temperature_on_leaf_growth_rate(16) == 0.2655638090621652
    assert temperature_on_leaf_growth_rate(20) == 0.5446031625628341
    assert temperature_on_leaf_growth_rate(24) == 0.7620184603088561
    assert temperature_on_leaf_growth_rate(27) == 0.9422215274064722
    assert temperature_on_leaf_growth_rate(30) == 1.0000351519210529
    assert temperature_on_leaf_growth_rate(36) == 0.7351624666288009
    assert temperature_on_leaf_growth_rate(42) == 0
