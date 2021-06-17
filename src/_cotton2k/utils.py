from datetime import date, datetime
from typing import Optional, Union


def date2doy(d: Union[str, int, date]) -> int:
    if isinstance(d, str):
        d = datetime.strptime(d, "%Y-%m-%d")
    result = 0
    if isinstance(d, date):
        result = d.timetuple().tm_yday
    elif isinstance(d, int) and d > 0:
        result = d
    return result


def doy2date(year: int, j: int) -> Optional[date]:
    try:
        return datetime.strptime(f"{year} {j}", "%Y %j").date()
    except ValueError:
        return None
