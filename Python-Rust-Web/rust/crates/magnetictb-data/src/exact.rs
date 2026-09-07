use serde::{Deserialize, Serialize};

use crate::DataError;

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ExactExpression {
    Integer {
        value: String,
    },
    Rational {
        numerator: String,
        denominator: String,
    },
    Symbol {
        name: String,
    },
    RootOfUnity {
        order: usize,
        power: i64,
    },
    Call {
        head: String,
        arguments: Vec<Self>,
    },
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct ExactMatrix {
    kind: MatrixTag,
    rows: usize,
    columns: usize,
    entries: Vec<ExactExpression>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
enum MatrixTag {
    #[serde(rename = "matrix")]
    Matrix,
}

impl ExactMatrix {
    #[must_use]
    pub const fn rows(&self) -> usize {
        self.rows
    }

    #[must_use]
    pub const fn columns(&self) -> usize {
        self.columns
    }

    #[must_use]
    pub fn entries(&self) -> &[ExactExpression] {
        &self.entries
    }

    pub(crate) fn validate(&self, field: &str) -> Result<(), DataError> {
        if self
            .rows
            .checked_mul(self.columns)
            .is_none_or(|length| length != self.entries.len())
        {
            return Err(DataError::new(
                "MalformedResource",
                format!("{field} dimensions do not match its row-major entries"),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(untagged)]
pub enum EncodedValue {
    Null(()),
    Boolean(bool),
    String(String),
    Integer(i64),
    Exact(ExactExpression),
    Matrix(ExactMatrix),
    List(EncodedList),
    OrderedAssociation(OrderedAssociation),
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct EncodedList {
    kind: ListTag,
    items: Vec<EncodedValue>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
enum ListTag {
    #[serde(rename = "list")]
    List,
}

impl EncodedList {
    #[must_use]
    pub fn items(&self) -> &[EncodedValue] {
        &self.items
    }
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct OrderedAssociation {
    kind: OrderedAssociationTag,
    entries: Vec<OrderedAssociationEntry>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
enum OrderedAssociationTag {
    #[serde(rename = "ordered_association")]
    OrderedAssociation,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct OrderedAssociationEntry {
    key: EncodedValue,
    value: EncodedValue,
}

impl OrderedAssociation {
    #[must_use]
    pub fn entries(&self) -> &[OrderedAssociationEntry] {
        &self.entries
    }
}

impl OrderedAssociationEntry {
    #[must_use]
    pub const fn key(&self) -> &EncodedValue {
        &self.key
    }

    #[must_use]
    pub const fn value(&self) -> &EncodedValue {
        &self.value
    }
}
