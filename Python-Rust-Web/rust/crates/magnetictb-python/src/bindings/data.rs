#[pyclass(name = "DataCatalog", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyDataCatalog {
    inner: Arc<CoreDataCatalog>,
}
#[pymethods]
impl PyDataCatalog {
    #[new]
    fn new() -> PyResult<Self> {
        CoreDataCatalog::load_embedded()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|error| python_data_error(&error))
    }

    fn counts_json(&self) -> PyResult<String> {
        serialize_data(&self.inner.counts())
    }

    fn manifest_json(&self) -> PyResult<String> {
        serialize_data(self.inner.manifest())
    }

    fn bravais_ids(&self) -> Vec<String> {
        self.inner
            .bravais_lattices()
            .iter()
            .map(|lattice| lattice.stable_id().to_owned())
            .collect()
    }

    fn msg_ids(&self) -> Vec<String> {
        self.inner
            .msg_groups()
            .iter()
            .map(|group| group.stable_id().to_owned())
            .collect()
    }

    fn rod_group_ids(&self) -> Vec<String> {
        self.inner
            .rod()
            .groups()
            .iter()
            .map(|group| group.stable_id().to_owned())
            .collect()
    }

    fn layer_group_ids(&self) -> Vec<String> {
        self.inner
            .layer()
            .groups()
            .iter()
            .map(|group| group.stable_id().to_owned())
            .collect()
    }

    fn bravais_json(&self, stable_id: &str) -> PyResult<String> {
        serialize_data(
            self.inner
                .bravais(stable_id)
                .ok_or_else(|| unknown_data_id("Bravais lattice", stable_id))?,
        )
    }

    fn msg_json(&self, stable_id: &str) -> PyResult<String> {
        serialize_data(
            self.inner
                .msg(stable_id)
                .ok_or_else(|| unknown_data_id("MSG", stable_id))?,
        )
    }

    fn resolve_msg_id_by_bns(&self, space_group: usize, serial: usize) -> PyResult<String> {
        self.inner
            .msg_by_bns(&[space_group, serial])
            .map(|group| group.stable_id().to_owned())
            .ok_or_else(|| {
                python_data_error(&CoreDataError::new(
                    "UnknownBnsKey",
                    format!("unknown BNS key ({space_group}, {serial})"),
                ))
            })
    }

    fn resolve_msg_id_by_source_index(&self, source_index: usize) -> PyResult<String> {
        self.inner
            .msg_by_source_index(source_index)
            .map(|group| group.stable_id().to_owned())
            .ok_or_else(|| {
                python_data_error(&CoreDataError::new(
                    "UnknownMagneticGroupIndex",
                    format!("unknown one-based magnetic-group source index {source_index}"),
                ))
            })
    }

    fn resolve_msg_id_by_classification_key(
        &self,
        classification: &str,
        source_key: Vec<usize>,
    ) -> PyResult<String> {
        let source_key = source_key.into_boxed_slice();
        self.inner
            .msg_by_classification_key(classification, &source_key)
            .map(|group| group.stable_id().to_owned())
            .ok_or_else(|| {
                python_data_error(&CoreDataError::new(
                    "UnknownMagneticGroupSelector",
                    format!("unknown {classification} selector key {source_key:?}"),
                ))
            })
    }

    fn resolve_msg_source_index_by_classification_key(
        &self,
        classification: &str,
        source_key: Vec<usize>,
    ) -> PyResult<usize> {
        let source_key = source_key.into_boxed_slice();
        self.inner
            .msg_by_classification_key(classification, &source_key)
            .map(CoreMsgGroup::source_index)
            .ok_or_else(|| {
                python_data_error(&CoreDataError::new(
                    "UnknownMagneticGroupSelector",
                    format!("unknown {classification} selector key {source_key:?}"),
                ))
            })
    }

    fn resolve_msg_id_by_classification(
        &self,
        classification: &str,
        source_index: usize,
    ) -> PyResult<String> {
        self.inner
            .msg_by_classification_index(classification, source_index)
            .map(|group| group.stable_id().to_owned())
            .ok_or_else(|| {
                python_data_error(&CoreDataError::new(
                    "UnknownMagneticGroupSelector",
                    format!("unknown {classification} selector at one-based index {source_index}"),
                ))
            })
    }

    fn wyckoff_json(&self, msg_id: &str) -> PyResult<String> {
        serialize_data(
            self.inner
                .wyckoff(msg_id)
                .ok_or_else(|| unknown_data_id("Wyckoff MSG", msg_id))?,
        )
    }

    fn resolve_wyckoff_json(
        &self,
        msg_id: &str,
        letter: &str,
        source_ordinal: Option<usize>,
    ) -> PyResult<String> {
        serialize_data(
            self.inner
                .wyckoff_entry(msg_id, letter, source_ordinal)
                .map_err(|error| python_data_error(&error))?,
        )
    }

    fn rod_group_json(&self, stable_id: &str) -> PyResult<String> {
        serialize_data(
            self.inner
                .rod_group(stable_id)
                .ok_or_else(|| unknown_data_id("rod group", stable_id))?,
        )
    }

    fn layer_group_json(&self, stable_id: &str) -> PyResult<String> {
        serialize_data(
            self.inner
                .layer_group(stable_id)
                .ok_or_else(|| unknown_data_id("layer group", stable_id))?,
        )
    }

    fn resolve_subperiodic_group_id_by_og(
        &self,
        kind: &str,
        og_key: Vec<usize>,
    ) -> PyResult<String> {
        let og_key = og_key.into_boxed_slice();
        self.inner
            .subperiodic_group_by_og(kind, &og_key)
            .map(|group| group.stable_id().to_owned())
            .ok_or_else(|| {
                python_data_error(&CoreDataError::new(
                    "UnknownSubperiodicGroupSelector",
                    format!("unknown {kind} group OG key {og_key:?}"),
                ))
            })
    }

    fn resolve_subperiodic_gray_key(&self, kind: &str, source_key: usize) -> PyResult<Vec<usize>> {
        self.inner
            .subperiodic_gray_key(kind, source_key)
            .ok_or_else(|| {
                python_data_error(&CoreDataError::new(
                    "UnknownSubperiodicGraySelector",
                    format!("unknown {kind} gray selector key {source_key}"),
                ))
            })
    }

    fn subperiodic_basic_vectors_json(&self, kind: &str, stable_id: &str) -> PyResult<String> {
        let group = match kind {
            "rod" => self.inner.rod_group(stable_id),
            "layer" => self.inner.layer_group(stable_id),
            _ => None,
        }
        .ok_or_else(|| unknown_data_id("subperiodic group", stable_id))?;
        serialize_data(
            self.inner
                .subperiodic_basic_vectors(kind, group.bravais_lattice_id())
                .ok_or_else(|| {
                    python_data_error(&CoreDataError::new(
                        "MissingSubperiodicBasicVectors",
                        format!(
                            "{} has no {} basic-vector record",
                            group.bravais_lattice_id(),
                            kind
                        ),
                    ))
                })?,
        )
    }

    fn classification_maps_json(&self) -> PyResult<String> {
        serialize_data(self.inner.classification_maps())
    }

    fn compile_msg_group_json(&self, stable_id: &str) -> PyResult<String> {
        let source = self
            .inner
            .msg(stable_id)
            .ok_or_else(|| unknown_data_id("MSG", stable_id))?;
        let compiled = compile_msg_group(source).map_err(|error| python_symmetry_error(&error))?;
        let context = compiled.operations()[0].spatial().rotation().context();
        let product_lattice_translations = compiled
            .product_records()
            .iter()
            .map(|row| {
                Value::Array(
                    row.iter()
                        .map(|record| {
                            Value::Array(
                                record
                                    .lattice_translation
                                    .iter()
                                    .map(element_to_json)
                                    .collect(),
                            )
                        })
                        .collect(),
                )
            })
            .collect::<Vec<_>>();
        serialize_json(&serde_json::json!({
            "stable_id": source.stable_id(),
            "classification": source.classification(),
            "context": context_to_json(context),
            "operations": compiled.operations().iter().map(|operation| serde_json::json!({
                "stable_id": operation.stable_id(),
                "label": operation.label(),
                "rotation": matrix_to_json(operation.spatial().rotation()),
                "translation": operation.spatial().translation().iter()
                    .map(element_to_json).collect::<Vec<_>>(),
                "antiunitary": operation.antiunitary()
            })).collect::<Vec<_>>(),
            "ordered_operation_ids": compiled.operations().iter().map(|operation| {
                operation.stable_id()
            }).collect::<Vec<_>>(),
            "ordered_operation_labels": compiled.operations().iter().map(|operation| {
                operation.label()
            }).collect::<Vec<_>>(),
            "antiunitary_flags": compiled.antiunitary_flags(),
            "multiplication_table": compiled.group().multiplication_table(),
            "identity": compiled.group().identity(),
            "inverse_indices": compiled.group().inverse_indices(),
            "product_lattice_translations": product_lattice_translations
        }))
        .map_err(|error| python_error(&error))
    }

    fn compile_msg_wyckoff_sites_json(&self, payload: &str) -> PyResult<String> {
        let request_value = parse_json(payload).map_err(|error| python_error(&error))?;
        let request = object(&request_value).map_err(|error| python_error(&error))?;
        let context = context_from_json(
            required(request, "context").map_err(|error| python_error(&error))?,
            128,
        )
        .map_err(|error| python_error(&error))?;
        let msg_id = required(request, "msg_id")
            .map_err(|error| python_error(&error))?
            .as_str()
            .ok_or_else(|| {
                python_error(&CoreExactError::new(
                    "MalformedSerialization",
                    "msg_id must be a string",
                ))
            })?;
        let source_ordinal = required(request, "source_ordinal")
            .map_err(|error| python_error(&error))?
            .as_u64()
            .and_then(|value| usize::try_from(value).ok())
            .ok_or_else(|| {
                python_error(&CoreExactError::new(
                    "MalformedSerialization",
                    "source_ordinal must be a nonnegative integer",
                ))
            })?;
        let encoded_bindings =
            object(required(request, "bindings").map_err(|error| python_error(&error))?)
                .map_err(|error| python_error(&error))?;
        let bindings = encoded_bindings
            .iter()
            .map(|(name, value)| {
                element_from_json(&context, value)
                    .map(|element| (name.clone(), element))
                    .map_err(|error| python_error(&error))
            })
            .collect::<PyResult<BTreeMap<_, _>>>()?;
        let compiled = compile_msg_wyckoff_sites(&self.inner, msg_id, source_ordinal, &bindings)
            .map_err(|error| python_error(&CoreExactError::new(error.tag(), error.to_string())))?;
        let sites = compiled.sites();
        let site_orbits = sites
            .site_orbits()
            .iter()
            .map(|orbit| {
                Value::Array(
                    orbit
                        .iter()
                        .map(|site| Value::Array(site.iter().map(element_to_json).collect()))
                        .collect(),
                )
            })
            .collect::<Vec<_>>();
        let cell_translations = sites
            .cell_translations()
            .iter()
            .map(|orbit| {
                Value::Array(
                    orbit
                        .iter()
                        .map(|operation| {
                            Value::Array(
                                operation
                                    .iter()
                                    .map(|translation| {
                                        Value::Array(
                                            translation.iter().map(element_to_json).collect(),
                                        )
                                    })
                                    .collect(),
                            )
                        })
                        .collect(),
                )
            })
            .collect::<Vec<_>>();
        serialize_json(&serde_json::json!({
            "msg_id": msg_id,
            "source_ordinal": compiled.source_ordinal(),
            "letter": compiled.letter(),
            "multiplication_table": compiled.symmetry().group().multiplication_table(),
            "antiunitary_flags": compiled.symmetry().antiunitary_flags(),
            "site_orbits": site_orbits,
            "image_site_indices": sites.image_site_indices(),
            "cell_translations": cell_translations
        }))
        .map_err(|error| python_error(&error))
    }

    fn compile_directed_bond_orbits_json(&self, payload: &str) -> PyResult<String> {
        geometry_boundary::compute_with_catalog(payload, &self.inner)
            .map_err(|error| python_error(&error))
    }

    fn compile_scalar_physical_representation_json(&self, payload: &str) -> PyResult<String> {
        physical_boundary::compute_with_catalog(payload, &self.inner)
            .map_err(|error| python_error(&error))
    }
}
