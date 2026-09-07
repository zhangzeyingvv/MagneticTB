use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExactError {
    tag: String,
    message: String,
}

impl ExactError {
    #[must_use]
    pub fn new(tag: impl Into<String>, message: impl Into<String>) -> Self {
        Self {
            tag: tag.into(),
            message: message.into(),
        }
    }

    #[must_use]
    pub fn tag(&self) -> &str {
        &self.tag
    }
}

impl Display for ExactError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}: {}", self.tag, self.message)
    }
}

impl Error for ExactError {}

pub type ExactResult<T> = Result<T, ExactError>;

pub(crate) fn fail<T>(tag: &str, message: impl Into<String>) -> ExactResult<T> {
    Err(ExactError::new(tag, message))
}
