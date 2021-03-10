from unittest.mock import patch

import pytest


def test_main():
    with patch("sys.argv", ["python", "-m", "cotton2k"]), pytest.raises(SystemExit):
        from cotton2k import __main__
