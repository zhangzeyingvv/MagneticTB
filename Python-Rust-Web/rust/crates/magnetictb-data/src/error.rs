use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DataError {
    tag: &'static str,
    message: String,
}

impl DataError {
    #[must_use]
    pub fn new(tag: &'static str, message: impl Into<String>) -> Self {
        Self {
            tag,
            message: message.into(),
        }
    }

    #[must_use]
    pub const fn tag(&self) -> &'static str {
        self.tag
    }
}

impl Display for DataError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}: {}", self.tag, self.message)
    }
}

impl Error for DataError {}
