import datetime
from typing import Optional


def strptime(date_string: str) -> datetime.date:
    return datetime.datetime.strptime(date_string, "%d-%b-%Y").date()


def date_to_day_of_year(
    date: datetime.date,
    start_year: Optional[int] = None,
) -> int:
    if start_year is None:
        start_year = date.year
    return (date - datetime.date(start_year, 1, 1)).days + 1
