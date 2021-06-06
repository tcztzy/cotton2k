from datetime import date, datetime
from typing import Union


def date2doy(d: Union[str, int, date]) -> int:
    if isinstance(d, str):
        d = datetime.strptime(d, "%Y-%m-%d")
    result = 0
    if isinstance(d, date):
        result = d.timetuple().tm_yday
    elif isinstance(d, int) and d > 0:
        result = d
    return result
