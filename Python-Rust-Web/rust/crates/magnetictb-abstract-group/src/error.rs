use std::error::Error;
use std::fmt::{self, Display, Formatter};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GroupError {
    tag: &'static str,
    stage: &'static str,
    detail: String,
}

impl GroupError {
    #[must_use]
    pub fn new(tag: &'static str, stage: &'static str, detail: impl Into<String>) -> Self {
        Self {
            tag,
            stage,
            detail: detail.into(),
        }
    }

    #[must_use]
    pub const fn tag(&self) -> &'static str {
        self.tag
    }

    #[must_use]
    pub const fn stage(&self) -> &'static str {
        self.stage
    }

    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl Display for GroupError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(formatter, "[{}:{}] {}", self.stage, self.tag, self.detail)
    }
}

impl Error for GroupError {}
