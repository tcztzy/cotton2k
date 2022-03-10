"""\
Cotton2K simulation model
=========================

```python
import cotton2k as c2k

c2k.run("/path/to/your/profile.toml")
```
"""
from os import PathLike

from .cotton2k import run as _run

__all__ = ("run",)


def run(profile: PathLike):
    _run(str(profile))
