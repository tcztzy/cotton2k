from datetime import date, datetime
from typing import Optional


def strptime(t) -> Optional[datetime.date]:
    t = t.strip()
    if len(t) == 11:
        return datetime.strptime(t, "%d-%b-%Y").date()


def date_to_day_of_year(d: date, start_year=None):
    if start_year is None:
        start_year = d.year
    return (d - date(start_year, 1, 1)).days + 1
