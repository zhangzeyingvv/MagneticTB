use std::error::Error;
use std::fmt::{Display, Formatter};

use cyclotomic_nullspace::ExactError;
use magnetictb_abstract_group::GroupError;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RepresentationError {
    tag: String,
    detail: String,
}

impl RepresentationError {
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
}

impl Display for RepresentationError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}: {}", self.tag, self.detail)
    }
}

impl Error for RepresentationError {}

impl From<ExactError> for RepresentationError {
    fn from(error: ExactError) -> Self {
        Self::new(error.tag(), error.to_string())
    }
}

impl From<GroupError> for RepresentationError {
    fn from(error: GroupError) -> Self {
        Self::new(error.tag(), error.to_string())
    }
}

pub type RepresentationResult<T> = Result<T, RepresentationError>;
