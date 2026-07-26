//! `#[derive(ScoreBlock)]` — projects a "score block" struct's fields into its
//! Parquet surface (the `ScoreBlock` trait's `column_schema`/`columns`), the
//! set-level ML feature-name walks (also on the trait), and the two ML feature
//! *value* lanes, which are INHERENT `[f64; N]` methods sized by inherent
//! `LINEAR_LEN`/`NONLINEAR_LEN` consts. All of it is walked off the one field
//! list, so the projections cannot desync the way a hand-maintained
//! `macro_rules!` invocation could.
//!
//! The lane values are inherent, not trait methods, on purpose: `fn f(&self) ->
//! [f64; Self::LEN]` in a trait needs unstable `generic_const_exprs`, while the
//! same signature on a concrete type is stable — and inherent consts compose
//! (`Fields::LINEAR_LEN = A::LINEAR_LEN + B::LINEAR_LEN`), which is what lets
//! the whole feature matrix be a compile-time width.
//!
//! See [`macro@ScoreBlock`] for the `#[feat(...)]` field grammar.

use proc_macro2::TokenStream;
use quote::{
    format_ident,
    quote,
};
use syn::spanned::Spanned;
use syn::{
    Data,
    DeriveInput,
    Error,
    Expr,
    Fields,
    Ident,
    Result,
    Type,
};

/// The scalar generators a `#[feat(...)]` field list may name.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Generator {
    Raw,
    Log2,
    Ln1p,
    Abs,
    Round,
    Isna,
}

impl Generator {
    fn parse(ident: &Ident) -> Result<Self> {
        match ident.to_string().as_str() {
            "raw" => Ok(Generator::Raw),
            "log2" => Ok(Generator::Log2),
            "ln1p" => Ok(Generator::Ln1p),
            "abs" => Ok(Generator::Abs),
            "round" => Ok(Generator::Round),
            "isna" => Ok(Generator::Isna),
            other => Err(Error::new(
                ident.span(),
                format!(
                    "unknown `#[feat(...)]` generator `{other}`; expected one of raw, log2, ln1p, abs, round, isna"
                ),
            )),
        }
    }

    /// Suffix appended to the bare field name to form the emitted feature name.
    ///
    /// This is the ONLY suffix table, and [`Generator::apply`] is the only
    /// arithmetic table; the two are read together at each generator site, so
    /// a name and the value under it are emitted from one match arm apiece.
    fn name_suffix(self) -> &'static str {
        match self {
            Generator::Raw => "",
            Generator::Log2 => "_log2",
            Generator::Ln1p => "_ln1p",
            Generator::Abs => "_abs",
            Generator::Round => "_round",
            Generator::Isna => "_isna",
        }
    }

    /// This generator's arithmetic, applied to an already-`f64` expression.
    /// `isna` is the only non-arithmetic one: a missingness indicator, so it is
    /// BY CONSTRUCTION never NaN, unlike `abs` which passes non-finites through.
    fn apply(self, x: TokenStream) -> TokenStream {
        match self {
            Generator::Raw => x,
            Generator::Log2 => quote! { (#x).log2() },
            Generator::Ln1p => quote! { (#x).ln_1p() },
            Generator::Abs => quote! { (#x).abs() },
            Generator::Round => quote! { (#x).round() },
            Generator::Isna => quote! { if (#x).is_finite() { 0.0 } else { 1.0 } },
        }
    }
}

/// Which feature lane a `#[feat(...)]` field's generators are routed into.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Lane {
    Linear,
    Nonlinear,
}

/// The scalar element types a field may declare.
#[derive(Clone, Copy)]
enum Scalar {
    F32,
    F64,
    U8,
    U32,
    Bool,
}

impl Scalar {
    fn from_type(ty: &Type) -> Option<Self> {
        let Type::Path(tp) = ty else { return None };
        let ident = tp.path.get_ident()?;
        match ident.to_string().as_str() {
            "f32" => Some(Scalar::F32),
            "f64" => Some(Scalar::F64),
            "u8" => Some(Scalar::U8),
            "u32" => Some(Scalar::U32),
            "bool" => Some(Scalar::Bool),
            _ => None,
        }
    }

    /// Method name shared by `SchemaSink`/`ColSink` for this scalar type
    /// (e.g. `SchemaSink::f32`, `ColSink::f32`).
    fn sink_method(self) -> Ident {
        let name = match self {
            Scalar::F32 => "f32",
            Scalar::F64 => "f64",
            Scalar::U8 => "u8",
            Scalar::U32 => "u32",
            Scalar::Bool => "bool",
        };
        format_ident!("{name}")
    }

    /// The scalar's own type, as tokens (e.g. `f32`) — used to reconstruct a
    /// field's type for the generated `sample()`'s `BlockFixture` turbofish.
    fn type_tokens(self) -> TokenStream {
        match self {
            Scalar::F32 => quote! { f32 },
            Scalar::F64 => quote! { f64 },
            Scalar::U8 => quote! { u8 },
            Scalar::U32 => quote! { u32 },
            Scalar::Bool => quote! { bool },
        }
    }

    /// `self.#ident` widened to the `f64` every feature-array slot holds.
    /// `bool` has no direct `as f64` cast (`error[E0606]`), so it goes through
    /// `u8` and lands as the usual 0.0/1.0 indicator.
    fn to_f64_expr(self, ident: &Ident) -> TokenStream {
        match self {
            Scalar::Bool => quote! { self.#ident as u8 as f64 },
            _ => quote! { self.#ident as f64 },
        }
    }
}

/// A field's shape: scalar, or a fixed-size `[f32; N]` array (`N` is an
/// arbitrary expression, emitted verbatim).
enum FieldShape {
    Scalar(Scalar),
    Array { len: Expr },
}

impl FieldShape {
    /// The field's type, as tokens — `f32`/`f64`/... for a scalar,
    /// `[f32; N]` for an array (reassembled from the parsed `len`).
    fn type_tokens(&self) -> TokenStream {
        match self {
            FieldShape::Scalar(scalar) => scalar.type_tokens(),
            FieldShape::Array { len } => quote! { [f32; #len] },
        }
    }

    fn from_type(ty: &Type) -> Result<Self> {
        if let Some(scalar) = Scalar::from_type(ty) {
            return Ok(FieldShape::Scalar(scalar));
        }
        if let Type::Array(arr) = ty {
            let Type::Path(elem) = arr.elem.as_ref() else {
                return Err(Error::new(
                    ty.span(),
                    "array fields must have element type `f32`",
                ));
            };
            if elem.path.get_ident().map(|i| i.to_string()).as_deref() != Some("f32") {
                return Err(Error::new(
                    ty.span(),
                    "array fields must have element type `f32`",
                ));
            }
            return Ok(FieldShape::Array {
                len: arr.len.clone(),
            });
        }
        Err(Error::new(
            ty.span(),
            "unsupported field type; expected one of f32, f64, u8, u32, bool, or [f32; N]",
        ))
    }
}

/// One generator as written in a `#[feat(...)]` list: the parsed variant plus
/// the token it came from, kept so diagnostics can echo the user's own
/// spelling and point at the offending token instead of the field name.
struct GeneratorSpec {
    generator: Generator,
    ident: Ident,
}

/// Parsed `#[feat(...)]` attribute contents: the requested generators plus
/// the optional `linear = <bool>` lane override.
struct FeatureAttr {
    generators: Vec<GeneratorSpec>,
    linear: bool,
}

/// Parses `#[feat(gen, gen, ..., linear = <bool>)]`: a comma-separated
/// list of distinct generator names plus an optional `linear = <bool>` item,
/// which may appear anywhere in the list.
fn parse_feature_attr(attr: &syn::Attribute) -> Result<FeatureAttr> {
    let mut generators: Vec<GeneratorSpec> = Vec::new();
    let mut linear: Option<bool> = None;

    attr.parse_nested_meta(|meta| {
        if meta.path.is_ident("linear") {
            if linear.is_some() {
                return Err(meta.error("duplicate `linear = ...` in #[feat(...)]"));
            }
            let value = meta.value()?;
            let lit: syn::LitBool = value.parse()?;
            linear = Some(lit.value);
            return Ok(());
        }
        let Some(ident) = meta.path.get_ident() else {
            return Err(meta.error("expected a generator name or `linear = <bool>`"));
        };
        let generator = Generator::parse(ident)?;
        // A repeat would emit the same feature name twice; the sinks would
        // happily accept the duplicate column.
        if generators.iter().any(|g| g.generator == generator) {
            return Err(Error::new(
                ident.span(),
                format!("duplicate `#[feat(...)]` generator `{ident}`"),
            ));
        }
        generators.push(GeneratorSpec {
            generator,
            ident: ident.clone(),
        });
        Ok(())
    })?;

    Ok(FeatureAttr {
        generators,
        linear: linear.unwrap_or(true),
    })
}

/// One named field of the derived struct, with its parsed shape and (if
/// present) its `#[feat(...)]` routing.
struct Field {
    ident: Ident,
    shape: FieldShape,
    /// `None` => column-only field (no `#[feat(...)]`).
    feature: Option<(Vec<GeneratorSpec>, Lane)>,
}

fn collect_fields(input: &DeriveInput) -> Result<Vec<Field>> {
    let Data::Struct(data) = &input.data else {
        return Err(Error::new(
            input.span(),
            "#[derive(ScoreBlock)] only supports structs with named fields",
        ));
    };
    let Fields::Named(named) = &data.fields else {
        return Err(Error::new(
            input.span(),
            "#[derive(ScoreBlock)] only supports structs with named fields",
        ));
    };

    let mut fields = Vec::with_capacity(named.named.len());
    for f in &named.named {
        let ident = f
            .ident
            .clone()
            .ok_or_else(|| Error::new(f.span(), "tuple struct fields are not supported"))?;
        let shape = FieldShape::from_type(&f.ty)?;

        let feature_attrs: Vec<&syn::Attribute> = f
            .attrs
            .iter()
            .filter(|a| a.path().is_ident("feat"))
            .collect();

        let feature = match feature_attrs.as_slice() {
            [] => None,
            [attr] => {
                let parsed = parse_feature_attr(attr)?;
                if parsed.generators.is_empty() {
                    return Err(Error::new(
                        attr.span(),
                        "#[feat(...)] requires at least one generator",
                    ));
                }
                let lane = if parsed.linear {
                    Lane::Linear
                } else {
                    Lane::Nonlinear
                };
                Some((parsed.generators, lane))
            }
            [_, dup, ..] => {
                return Err(Error::new(
                    dup.span(),
                    "duplicate #[feat(...)] attribute on field",
                ));
            }
        };

        fields.push(Field {
            ident,
            shape,
            feature,
        });
    }
    Ok(fields)
}

/// Emits the `SchemaSink`/`ColSink` calls shared by `column_schema` and
/// `columns` (every field is a parquet column, feature-routed or not).
fn schema_and_column_calls(fields: &[Field]) -> (Vec<TokenStream>, Vec<TokenStream>) {
    let mut schema_calls = Vec::with_capacity(fields.len());
    let mut column_calls = Vec::with_capacity(fields.len());

    for field in fields {
        let ident = &field.ident;
        let name = ident.to_string();
        match &field.shape {
            FieldShape::Scalar(scalar) => {
                let method = scalar.sink_method();
                schema_calls.push(quote! { out.#method(#name); });
                column_calls.push(quote! { out.#method(#name, self.#ident); });
            }
            FieldShape::Array { len } => {
                schema_calls.push(quote! { out.f32_array(#name, #len); });
                column_calls.push(quote! { out.f32_array(#name, &self.#ident); });
            }
        }
    }

    (schema_calls, column_calls)
}

/// One feature lane's codegen: the width const's terms, the statements that
/// fill the `[f64; LEN]` value array, and the set-level name pushes.
struct LaneCode {
    /// Summands of the lane's `LEN` const — `1usize` per scalar generator,
    /// the array's own length expression (verbatim, so a const path like
    /// `NUM_MS2_IONS` stays a const path) per array generator.
    len_terms: Vec<TokenStream>,
    /// Writes into `out` at the running `at` cursor; each advances `at` by
    /// exactly the term it contributed to `len_terms`, in the same order.
    fills: Vec<TokenStream>,
    /// `NameSink` pushes, in the same order as `fills`.
    name_calls: Vec<TokenStream>,
}

/// Emits the value-array fills, the width terms, and the `NameSink` calls for
/// one feature lane. All three walk the field list once, in one loop, so they
/// are index-for-index aligned by construction.
fn lane_code(fields: &[Field], lane: Lane) -> Result<LaneCode> {
    let mut code = LaneCode {
        len_terms: Vec::new(),
        fills: Vec::new(),
        name_calls: Vec::new(),
    };

    for field in fields {
        let Some((generators, field_lane)) = &field.feature else {
            continue;
        };
        if *field_lane != lane {
            continue;
        }
        let ident = &field.ident;
        let name = ident.to_string();

        for spec in generators {
            let generator = spec.generator;
            match &field.shape {
                FieldShape::Scalar(scalar) => {
                    let value = generator.apply(scalar.to_f64_expr(ident));
                    // One name, computed once here, used by the name walk; the
                    // value walk is positional, so it cannot name anything
                    // differently — there is nothing to name.
                    let feat_name = format!("{name}{}", generator.name_suffix());
                    code.len_terms.push(quote! { 1usize });
                    code.fills.push(quote! {
                        out[at] = #value;
                        at += 1;
                    });
                    code.name_calls.push(quote! { out.push(#feat_name); });
                }
                FieldShape::Array { len } => {
                    let (elem, suffix) = match generator {
                        Generator::Raw => (quote! { *v as f64 }, ""),
                        Generator::Isna => {
                            (quote! { if v.is_finite() { 0.0 } else { 1.0 } }, "_isna")
                        }
                        _ => {
                            let written = &spec.ident;
                            return Err(Error::new(
                                written.span(),
                                format!(
                                    "generator `{written}` is not supported on array fields; only raw and isna apply"
                                ),
                            ));
                        }
                    };
                    // Field names are identifiers, so the only brace in this
                    // format string is the `{i}` the generated loop fills in.
                    let name_fmt = format!("{name}_{{i}}{suffix}");
                    code.len_terms.push(quote! { (#len) });
                    code.fills.push(quote! {
                        for (k, v) in self.#ident.iter().enumerate() {
                            out[at + k] = #elem;
                        }
                        at += #len;
                    });
                    code.name_calls.push(quote! {
                        for i in 0..#len {
                            out.push(&format!(#name_fmt));
                        }
                    });
                }
            }
        }
    }

    Ok(code)
}

/// The inherent `pub const #len_const: usize` + `pub fn #method(&self) ->
/// [f64; Self::#len_const]` pair for one lane. An empty lane skips the cursor
/// entirely (a `let mut at` no statement ever touches is an `unused_mut`).
fn lane_methods(code: &LaneCode, len_const: Ident, method: Ident) -> TokenStream {
    let LaneCode {
        len_terms, fills, ..
    } = code;
    let body = if fills.is_empty() {
        quote! { [0.0f64; Self::#len_const] }
    } else {
        quote! {
            let mut out = [0.0f64; Self::#len_const];
            let mut at = 0usize;
            #(#fills)*
            debug_assert_eq!(at, Self::#len_const);
            out
        }
    };
    quote! {
        /// This lane's feature count — the sum of one per scalar generator and
        /// `N` per `[f32; N]` generator, so the lane's width is a compile-time
        /// constant that composes into the whole matrix's width.
        pub const #len_const: usize = 0usize #( + #len_terms )*;

        /// This lane's feature *values*, positionally aligned with the
        /// corresponding `*_feature_names` walk.
        pub fn #method(&self) -> [f64; Self::#len_const] {
            #body
        }
    }
}

/// Emits one `#field: <#ty as BlockFixture>::fixture()` initializer per
/// field, for the generated `sample()`.
fn sample_field_calls(fields: &[Field]) -> Vec<TokenStream> {
    fields
        .iter()
        .map(|field| {
            let ident = &field.ident;
            let ty = field.shape.type_tokens();
            quote! {
                #ident: <#ty as crate::scoring::blocks::BlockFixture>::fixture()
            }
        })
        .collect()
}

/// Pure codegen entry point: parses the `#[feat(...)]` grammar off each
/// named field of `input` and emits an `impl ScoreBlock for $Name` (Parquet
/// schema + columns, and the two set-level name walks) plus an inherent impl
/// carrying the two lane width consts, the two `[f64; N]` lane value methods,
/// and `sample()`. `sample()` stays un-gated (not `#[cfg(test)]`): it is called
/// uniformly across every block by `ScoringFields::sample_default`, which
/// itself must stay reachable in a normal (non-test) build because
/// `ScoringFields::sample` is used from another crate's test module
/// (`timsseek_cli`), compiled against timsseek's normal build. Unit-testable
/// without a proc-macro context — see `tests` below.
pub(crate) fn derive_score_block(input: DeriveInput) -> Result<TokenStream> {
    let fields = collect_fields(&input)?;
    let name = &input.ident;
    let (generics_impl, generics_ty, generics_where) = input.generics.split_for_impl();

    let (schema_calls, column_calls) = schema_and_column_calls(&fields);
    let linear = lane_code(&fields, Lane::Linear)?;
    let nonlinear = lane_code(&fields, Lane::Nonlinear)?;
    let linear_name_calls = &linear.name_calls;
    let nonlinear_name_calls = &nonlinear.name_calls;
    let linear_methods = lane_methods(
        &linear,
        format_ident!("LINEAR_LEN"),
        format_ident!("linear_feature_array"),
    );
    let nonlinear_methods = lane_methods(
        &nonlinear,
        format_ident!("NONLINEAR_LEN"),
        format_ident!("nonlinear_feature_array"),
    );
    let sample_calls = sample_field_calls(&fields);

    Ok(quote! {
        impl #generics_impl crate::scoring::blocks::ScoreBlock for #name #generics_ty #generics_where {
            fn column_schema(out: &mut crate::scoring::blocks::SchemaSink) {
                #(#schema_calls)*
            }

            fn columns(&self, out: &mut crate::scoring::blocks::ColSink) {
                #(#column_calls)*
            }

            fn linear_feature_names(out: &mut crate::scoring::blocks::NameSink) {
                #(#linear_name_calls)*
            }

            fn nonlinear_feature_names(out: &mut crate::scoring::blocks::NameSink) {
                #(#nonlinear_name_calls)*
            }
        }

        impl #generics_impl #name #generics_ty #generics_where {
            #linear_methods
            #nonlinear_methods

            /// Fixture with every field a fixed constant (see
            /// [`crate::scoring::blocks::BlockFixture`]). Transition
            /// scaffolding: un-gated because `ScoringFields::sample_default`
            /// is reachable from production code today; a later task moves
            /// this behind `cfg(test)` once that dependency is gone.
            pub fn sample() -> Self {
                Self {
                    #(#sample_calls),*
                }
            }
        }
    })
}

/// Derives `timsseek`'s `ScoreBlock` trait for a struct of named fields, so a
/// score family's Parquet projection and its ML feature projections are all
/// generated from the one field list that defines the struct.
///
/// # Fields and names
///
/// A field's **name is its name**: field `apex_lazyscore` becomes Parquet
/// column `apex_lazyscore` and (if featurized) feature `apex_lazyscore`, with
/// a per-generator suffix appended. There is no rename attribute.
///
/// Every field becomes a Parquet column. Fields *without* `#[feat(...)]` are
/// Parquet-column-only: they are emitted by `column_schema`/`columns` and
/// never reach a feature lane. `#[feat(...)]` is what additionally routes a
/// field into the ML matrix.
///
/// Supported field types: `f32`, `f64`, `u8`, `u32`, `bool`, and `[f32; N]`
/// (`N` may be any const expression; it is emitted verbatim).
///
/// # `#[feat(...)]` grammar
///
/// ```text
/// #[feat(<generator>, <generator>, ..., linear = <bool>)]
/// ```
///
/// At least one generator is required; repeats are rejected (they would emit
/// the same feature name twice), as is a second `#[feat(...)]` on one field.
/// The generators, and the suffix each appends to the field name:
///
/// | generator | feature name      | value                                  |
/// |-----------|-------------------|----------------------------------------|
/// | `raw`     | `{field}`         | the value itself                       |
/// | `log2`    | `{field}_log2`    | `x.log2()`                             |
/// | `ln1p`    | `{field}_ln1p`    | `x.ln_1p()`                            |
/// | `abs`     | `{field}_abs`     | `x.abs()` (NaN preserved)              |
/// | `round`   | `{field}_round`   | `x.round()`                            |
/// | `isna`    | `{field}_isna`    | `1.0` if non-finite, else `0.0`        |
///
/// Values are pushed as `f64`. A `bool` field is widened through `u8`, so it
/// lands as the usual 0.0/1.0 indicator.
///
/// ## Lanes
///
/// `linear = true` is the **default**: a plain `#[feat(raw)]` field goes to the
/// LINEAR lane (the monotone-transform lane LDA reads). The LDA lane is
/// therefore **opt-out**, not opt-in — write `linear = false` to route a field
/// to the NONLINEAR lane instead (the tree/GBM lane, which needs no monotone
/// transforms). A field is in exactly one lane; all of its generators go there
/// together.
///
/// ## Array fields
///
/// `[f32; N]` fields accept only `raw` and `isna` — any other generator is a
/// compile error pointing at the generator token. They fan out to one feature
/// per element, `{field}_0 .. {field}_{N-1}` (and `{field}_{i}_isna` for
/// `isna`), matching the `{field}_{i}` naming their Parquet columns already
/// use.
///
/// # What is generated
///
/// Four `ScoreBlock` trait methods:
///
/// - `column_schema(&mut SchemaSink)` — Parquet field dtypes/nullability
/// - `columns(&self, &mut ColSink)` — Parquet values
/// - `linear_feature_names(&mut NameSink)` / `nonlinear_feature_names(&mut NameSink)`
///
/// and, on an INHERENT impl, the lane value surface:
///
/// - `pub const LINEAR_LEN: usize` / `pub const NONLINEAR_LEN: usize`
/// - `pub fn linear_feature_array(&self) -> [f64; Self::LINEAR_LEN]`
/// - `pub fn nonlinear_feature_array(&self) -> [f64; Self::NONLINEAR_LEN]`
///
/// The lane values are inherent because `-> [f64; Self::LEN]` in a *trait*
/// needs unstable `generic_const_exprs`. Inherent, the width is a plain
/// compile-time constant that composes: a set of blocks sums their consts to
/// get the whole feature matrix's width, with no runtime length bookkeeping.
///
/// The name walk stays dynamic (`NameSink`) — an array field's `{field}_{i}`
/// fan-out is a loop over a const *path*, which no macro can unroll.
///
/// A lane's value array and its name walk are emitted from the same single
/// pass over the field list, so entry `j` of the array is always the feature
/// named `j` — they cannot desync.
///
/// Plus a public inherent `fn sample() -> Self` on the struct, filling every
/// field from `BlockFixture::fixture()` — a fixed, finite, non-zero constant
/// per type. It is deliberately not `#[cfg(test)]`; see `derive_score_block`.
///
/// # Example
///
/// ```ignore
/// #[derive(Debug, Clone, Copy, ScoreBlock)]
/// pub struct MyScores {
///     /// LINEAR lane (default), feature `intensity_log2`.
///     #[feat(log2)]
///     pub intensity: f32,
///     /// NONLINEAR lane, features `charge` and `charge_isna`.
///     #[feat(raw, isna, linear = false)]
///     pub charge: f32,
///     /// LINEAR lane, features `ion_err_0 .. ion_err_5`.
///     #[feat(raw)]
///     pub ion_err: [f32; 6],
///     /// Parquet column only — no feature.
///     pub n_cycles: u32,
/// }
/// ```
///
/// gives `MyScores::LINEAR_LEN == 1usize + (6)` and `MyScores::NONLINEAR_LEN ==
/// 1usize + 1usize`.
#[proc_macro_derive(ScoreBlock, attributes(feat))]
pub fn score_block_derive(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let di = syn::parse_macro_input!(input as DeriveInput);
    derive_score_block(di)
        .unwrap_or_else(|e| e.to_compile_error())
        .into()
}

#[cfg(test)]
mod tests {
    use super::derive_score_block;
    use quote::quote;

    #[test]
    fn emits_the_trait_and_inherent_lane_surface() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T {
                #[feat(log2)] pub a: f32,
                #[feat(raw, linear = false)] pub b: f32,
                pub c: f32,
            }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        assert!(ts.contains("fn column_schema"));
        assert!(ts.contains("fn columns"));
        assert!(ts.contains("fn linear_feature_names"));
        assert!(ts.contains("fn nonlinear_feature_names"));
        // The lane VALUES are inherent const-sized arrays, not trait methods.
        assert!(!ts.contains("fn linear_features"), "{ts}");
        assert!(!ts.contains("fn nonlinear_features"), "{ts}");
        let f = flat(&ts);
        assert!(
            f.contains("pubconstLINEAR_LEN:usize=0usize+1usize;"),
            "{ts}"
        );
        assert!(
            f.contains("pubconstNONLINEAR_LEN:usize=0usize+1usize;"),
            "{ts}"
        );
        assert!(
            f.contains("pubfnlinear_feature_array(&self)->[f64;Self::LINEAR_LEN]"),
            "{ts}"
        );
        // a (log2) -> linear lane; b (raw) -> nonlinear lane; c -> column-only
        assert!(ts.contains("log2"));
        // all three fields appear in column_schema
        assert!(ts.contains("\"c\""));
    }

    /// `TokenStream::to_string` spaces tokens out; the emitted *shape* is what
    /// matters, so compare whitespace-stripped.
    fn flat(ts: &str) -> String {
        ts.chars().filter(|c| !c.is_whitespace()).collect()
    }

    /// A lane with no fields must still emit its (zero) const and a method
    /// returning `[f64; 0]`, and must NOT open a cursor no statement touches
    /// (that is an `unused_mut` at the expansion site).
    #[test]
    fn empty_lane_emits_a_zero_length_array_without_a_cursor() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(raw)] pub a: f32 }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        let f = flat(&ts);
        assert!(f.contains("pubconstNONLINEAR_LEN:usize=0usize;"), "{ts}");
        assert!(
            f.contains(
                "pubfnnonlinear_feature_array(&self)->[f64;Self::NONLINEAR_LEN]\
                 {[0.0f64;Self::NONLINEAR_LEN]}"
            ),
            "{ts}"
        );
    }

    /// The name walk owns the only copy of a feature's name; the value walk is
    /// positional. `Generator::name_suffix` is still the only suffix table, and
    /// `Generator::apply` the only arithmetic table.
    #[test]
    fn scalar_name_is_suffixed_and_the_value_is_positional() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(log2)] pub a: f32 }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        let f = flat(&ts);
        assert!(f.contains(r#"out[at]=(self.aasf64).log2();at+=1;"#), "{ts}");
        assert!(f.contains(r#"out.push("a_log2");"#), "{ts}");
        // The bare `"a"` survives only as the Parquet column (schema + value);
        // `"a_log2"` exists exactly once, in the name walk.
        assert_eq!(f.matches(r#""a""#).count(), 2, "{ts}");
        assert_eq!(f.matches(r#""a_log2""#).count(), 1, "{ts}");
    }

    /// `Generator::Raw`'s empty suffix must still yield the bare name, and a
    /// multi-generator field must emit one slot per generator, in written
    /// order, with the array width counting both.
    #[test]
    fn raw_keeps_bare_name_and_generators_keep_order() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(raw, isna)] pub x: f32 }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        let f = flat(&ts);
        assert!(
            f.contains(
                "out[at]=self.xasf64;at+=1;\
                 out[at]=if(self.xasf64).is_finite(){0.0}else{1.0};at+=1;"
            ),
            "{ts}"
        );
        assert!(f.contains(r#"out.push("x");out.push("x_isna");"#), "{ts}");
        assert!(
            f.contains("pubconstLINEAR_LEN:usize=0usize+1usize+1usize;"),
            "{ts}"
        );
    }

    #[test]
    fn emits_sample() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T {
                #[feat(raw)] pub a: f32,
                pub arr: [f32; 3],
            }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        assert!(ts.contains("fn sample"));
        assert!(ts.contains("BlockFixture"));
        assert!(ts.contains("fixture"));
    }

    #[test]
    fn rejects_unknown_generator() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(bogus)] pub a: f32 }
        })
        .unwrap();
        assert!(derive_score_block(di).is_err());
    }

    /// An array field's length is a const PATH, so it must reach the width
    /// const verbatim (never a resolved literal, which the macro cannot know)
    /// and the fill must advance the cursor by that same path.
    #[test]
    fn array_field_width_is_the_const_path_verbatim() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(raw)] pub arr: [f32; NUM_MS2_IONS] }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        let f = flat(&ts);
        assert!(ts.contains("f32_array"));
        assert!(
            f.contains("pubconstLINEAR_LEN:usize=0usize+(NUM_MS2_IONS);"),
            "{ts}"
        );
        assert!(
            f.contains(
                "for(k,v)inself.arr.iter().enumerate(){out[at+k]=*vasf64;}at+=NUM_MS2_IONS;"
            ),
            "{ts}"
        );
        // The `{field}_{i}` name fan-out stays a runtime loop — no macro can
        // unroll a const path — but the format string itself is a literal.
        assert!(
            f.contains(r#"foriin0..NUM_MS2_IONS{out.push(&format!("arr_{i}"));}"#),
            "{ts}"
        );
    }

    /// The array `isna` companion suffixes the INDEX, not the field: the
    /// columns are `{field}_{i}_isna`, matching what the Parquet fan-out and
    /// the pre-array-derive name table both used.
    #[test]
    fn array_isna_names_suffix_the_index() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(isna)] pub arr: [f32; NUM_MS2_IONS] }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        assert!(
            flat(&ts).contains(r#"out.push(&format!("arr_{i}_isna"));"#),
            "{ts}"
        );
    }

    #[test]
    fn rejects_non_raw_generator_on_array() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(log2)] pub arr: [f32; NUM_MS2_IONS] }
        })
        .unwrap();
        let err = derive_score_block(di).unwrap_err();
        // Echoes the spelling the user wrote, not the `Generator` variant.
        let msg = err.to_string();
        assert!(msg.contains("`log2`"), "{msg}");
        assert!(!msg.contains("Log2"), "{msg}");
    }

    /// `bool as f64` is not a legal cast (E0606); the widening must go through
    /// `u8` or the expansion fails to compile at the *use* site.
    #[test]
    fn bool_field_widens_through_u8() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(raw)] pub flag: bool }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        assert!(ts.contains("as u8 as f64"), "{ts}");
    }

    #[test]
    fn rejects_duplicate_generator() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(raw, raw)] pub a: f32 }
        })
        .unwrap();
        let err = derive_score_block(di).unwrap_err();
        assert!(err.to_string().contains("duplicate"), "{err}");
    }
}
