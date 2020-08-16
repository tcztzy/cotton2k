from datetime import date, datetime


def strptime(t) -> date:
    return datetime.strptime(t, "%d-%b-%Y").date()


def date_to_day_of_year(d: date, start_year=None) -> int:
    if start_year is None:
        start_year = d.year
    return (d - date(start_year, 1, 1)).days + 1
