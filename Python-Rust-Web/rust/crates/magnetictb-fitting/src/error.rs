use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FittingError {
    tag: String,
    detail: String,
}

impl FittingError {
    #[must_use]
    pub fn new(tag: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            tag: tag.into(),
            detail: detail.into(),
        }
    }

    #[must_use]
    pub fn tag(&self) -> &str {
        &self.tag
    }

    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl Display for FittingError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}: {}", self.tag, self.detail)
    }
}

impl Error for FittingError {}

pub type FittingResult<T> = Result<T, FittingError>;
