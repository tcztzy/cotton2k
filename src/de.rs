struct NaiveDateVisitor;

impl<'de> serde::de::Visitor<'de> for NaiveDateVisitor {
    type Value = chrono::NaiveDate;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "a string represents chrono::NaiveDate")
    }

    fn visit_str<E>(self, s: &str) -> Result<Self::Value, E>
    where
        E: serde::de::Error,
    {
        match chrono::NaiveDate::parse_from_str(s, "%F") {
            Ok(t) => Ok(t),
            Err(_) => Err(serde::de::Error::invalid_value(serde::de::Unexpected::Str(s), &self)),
        }
    }
}
struct OptionalNaiveDateVisitor;

impl<'de> serde::de::Visitor<'de> for OptionalNaiveDateVisitor {
    type Value = Option<chrono::NaiveDate>;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "a optional string represents chrono::NaiveDate")
    }

    fn visit_none<E>(self) -> Result<Self::Value, E>
    where
        E: serde::de::Error,
    {
        Ok(None)
    }

    fn visit_some<D>(self, deserializer: D) -> Result<Self::Value, D::Error>
    where
        D: serde::de::Deserializer<'de>,
    {
        Ok(Some(deserializer.deserialize_str(NaiveDateVisitor)?))
    }
}

pub fn from_isoformat<'de, D>(d: D) -> Result<chrono::NaiveDate, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    d.deserialize_str(NaiveDateVisitor)
}

pub fn from_isoformat_option<'de, D>(d: D) -> Result<Option<chrono::NaiveDate>, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    d.deserialize_option(OptionalNaiveDateVisitor)
}
