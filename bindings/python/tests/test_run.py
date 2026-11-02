"""Exercise the public Python API against the installed native extension."""

import csv
import json
import shutil
from pathlib import Path

import pytest

import cotton2k

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


@pytest.fixture
def profile(tmp_path):
    fixture_dir = Path(__file__).resolve().parents[3] / "tests/fixtures/minimal"
    scenario = tmp_path / "scenario"
    scenario.mkdir()
    for name in ("profile.toml", "weather.csv", "soil_imp.csv"):
        shutil.copyfile(fixture_dir / name, scenario / name)
    path = scenario / "profile.toml"
    path.write_text(
        path.read_text().replace('stop_date = "2020-03-31"', 'stop_date = "2020-04-03"')
    )
    return path


@pytest.mark.parametrize("path_type", [str, Path], ids=["str", "pathlib"])
@pytest.mark.parametrize("extension", ["toml", "json"])
def test_run_three_day_simulation(profile, tmp_path, monkeypatch, path_type, extension):
    if extension == "json":
        contents = tomllib.loads(profile.read_text())
        profile = profile.with_suffix(".json")
        profile.write_text(json.dumps(contents))

    # Input CSV paths must resolve beside the profile, independently of cwd.
    monkeypatch.chdir(tmp_path)
    assert cotton2k.run(path_type(profile.relative_to(tmp_path))) is None

    with (profile.parent / "output.csv").open(newline="") as output:
        rows = list(csv.DictReader(output))
    assert [row["date"] for row in rows] == ["2020-04-01", "2020-04-02", "2020-04-03"]

    # Reference values are shared with the Rust three-day regression test.
    expected_water = [
        (0.37218740888258084, 0.36759272052047426, 0.3665730085273955),
        (0.3705127990396385, 0.36786902980003505, 0.3662366712173369),
        (0.3536694472153136, 0.3677432905162993, 0.3660897425976588),
    ]
    for row, water in zip(rows, expected_water):
        assert float(row["leaf_area_index"]) == pytest.approx(0.001, rel=0, abs=1e-10)
        assert float(row["root_weight"]) == pytest.approx(0.208, rel=0, abs=1e-10)
        assert [float(row[column]) for column in ("swc0-10", "swc0-20", "swc0-30")] == (
            pytest.approx(water, rel=0, abs=1e-10)
        )

    metadata = json.loads((profile.parent / "meta.json").read_text())
    assert metadata["status"] == "succeeded"
    assert metadata["error"] is None
    assert (profile.parent / "input" / profile.name).read_text() == profile.read_text()
    assert not (tmp_path / "output.csv").exists()


@pytest.mark.parametrize("path_type", [str, Path], ids=["str", "pathlib"])
def test_missing_profile_raises_value_error(tmp_path, path_type):
    with pytest.raises(ValueError, match="profile file.*does not exist"):
        cotton2k.run(path_type(tmp_path / "missing.toml"))


@pytest.mark.parametrize(
    ("extension", "contents", "message"),
    [
        ("toml", "latitude = [", "failed to parse TOML profile"),
        ("json", "{", "failed to parse JSON profile"),
        ("txt", "", "unsupported profile extension"),
    ],
)
def test_invalid_profile_raises_value_error(tmp_path, extension, contents, message):
    path = tmp_path / f"profile.{extension}"
    path.write_text(contents)
    with pytest.raises(ValueError, match=message):
        cotton2k.run(path)


def test_invalid_date_range_raises_value_error(profile):
    profile.write_text(
        profile.read_text().replace('stop_date = "2020-04-03"', 'stop_date = "2020-03-31"')
    )
    with pytest.raises(ValueError, match="stop_date.*must not be earlier than start_date"):
        cotton2k.run(profile)


@pytest.mark.parametrize("name", ["weather.csv", "soil_imp.csv"])
def test_missing_input_raises_value_error(profile, name):
    (profile.parent / name).unlink()
    with pytest.raises(ValueError, match="does not exist"):
        cotton2k.run(profile)


@pytest.mark.parametrize(
    ("original", "invalid", "message"),
    [
        ("2020-04-02,19.0", "2020-04-02,not-a-number", "failed to parse weather file"),
        ("2020-04-02", "2022-04-02", "weather record date 2022-04-02"),
    ],
)
def test_invalid_weather_raises_value_error(profile, original, invalid, message):
    weather = profile.parent / "weather.csv"
    weather.write_text(weather.read_text().replace(original, invalid))
    with pytest.raises(ValueError, match=message):
        cotton2k.run(profile)


def test_artifact_io_failure_raises_runtime_error(profile):
    (profile.parent / "input").write_text("blocks the snapshot directory")
    with pytest.raises(RuntimeError, match="failed to create input snapshot directory"):
        cotton2k.run(profile)
