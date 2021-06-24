from unittest.mock import patch

import pytest


def test_main():
    with patch("sys.argv", ["cotton2k"]), pytest.raises(SystemExit):
        from cotton2k import __main__


def test_file_not_found():
    with patch("sys.argv", ["cotton2k", "NOT_EXIST"]), pytest.raises(FileNotFoundError):
        from cotton2k import __main__
