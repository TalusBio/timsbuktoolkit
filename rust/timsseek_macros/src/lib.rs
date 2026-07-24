//! `#[derive(ScoreBlock)]` — projects a "score block" struct's fields into the
//! six `ScoreBlock` trait methods (schema, columns, and the two ML feature
//! lanes), so the projections cannot desync the way a hand-maintained
//! `score_block!` macro_rules invocation could.
//!
//! See the field grammar in [`derive_score_block`].

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

    /// `FrameSink` method used to push a scalar value for this generator.
    fn scalar_push_method(self) -> Ident {
        let name = match self {
            Generator::Raw => "push",
            Generator::Log2 => "push_log2",
            Generator::Ln1p => "push_ln1p",
            Generator::Abs => "push_abs",
            Generator::Round => "push_round",
            Generator::Isna => "push_isna",
        };
        format_ident!("{name}")
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
}

/// A field's shape: scalar, or a fixed-size `[f32; N]` array (`N` is an
/// arbitrary expression, emitted verbatim).
enum FieldShape {
    Scalar(Scalar),
    Array { len: Expr },
}

impl FieldShape {
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

/// Parsed `#[feat(...)]` attribute contents: the requested generators plus
/// the optional `linear = <bool>` lane override.
struct FeatureAttr {
    generators: Vec<Generator>,
    linear: bool,
}

/// Parses `#[feat(gen, gen, ..., linear = <bool>)]`: a comma-separated
/// list of generator names plus an optional `linear = <bool>` item, which may
/// appear anywhere in the list.
fn parse_feature_attr(attr: &syn::Attribute) -> Result<FeatureAttr> {
    let mut generators = Vec::new();
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
        generators.push(Generator::parse(ident)?);
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
    feature: Option<(Vec<Generator>, Lane)>,
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

/// Emits the `FrameSink`/`NameSink` calls for one feature lane.
fn lane_calls(fields: &[Field], lane: Lane) -> Result<(Vec<TokenStream>, Vec<TokenStream>)> {
    let mut feature_calls = Vec::new();
    let mut name_calls = Vec::new();

    for field in fields {
        let Some((generators, field_lane)) = &field.feature else {
            continue;
        };
        if *field_lane != lane {
            continue;
        }
        let ident = &field.ident;
        let name = ident.to_string();

        for &generator in generators {
            match &field.shape {
                FieldShape::Scalar(_) => {
                    let push = generator.scalar_push_method();
                    feature_calls.push(quote! {
                        out.#push(#name, self.#ident as f64);
                    });
                    let feat_name = format!("{name}{}", generator.name_suffix());
                    name_calls.push(quote! { out.push(#feat_name); });
                }
                FieldShape::Array { len } => match generator {
                    Generator::Raw => {
                        feature_calls.push(quote! {
                            out.push_slice(#name, &self.#ident);
                        });
                        name_calls.push(quote! {
                            for i in 0..#len {
                                out.push(&format!("{}_{i}", #name));
                            }
                        });
                    }
                    Generator::Isna => {
                        feature_calls.push(quote! {
                            out.push_slice_isna(#name, &self.#ident);
                        });
                        name_calls.push(quote! {
                            for i in 0..#len {
                                out.push(&format!("{}_{i}_isna", #name));
                            }
                        });
                    }
                    other => {
                        return Err(Error::new(
                            ident.span(),
                            format!(
                                "generator `{other:?}` is not supported on array fields; only raw and isna apply"
                            ),
                        ));
                    }
                },
            }
        }
    }

    Ok((feature_calls, name_calls))
}

/// Pure codegen entry point: parses the `#[feat(...)]` grammar off each
/// named field of `input` and emits an `impl ScoreBlock for $Name` with the
/// six projection methods described in the crate docs. Unit-testable without
/// a proc-macro context — see `tests` below.
pub(crate) fn derive_score_block(input: DeriveInput) -> Result<TokenStream> {
    let fields = collect_fields(&input)?;
    let name = &input.ident;
    let (generics_impl, generics_ty, generics_where) = input.generics.split_for_impl();

    let (schema_calls, column_calls) = schema_and_column_calls(&fields);
    let (linear_feature_calls, linear_name_calls) = lane_calls(&fields, Lane::Linear)?;
    let (nonlinear_feature_calls, nonlinear_name_calls) = lane_calls(&fields, Lane::Nonlinear)?;

    Ok(quote! {
        impl #generics_impl crate::scoring::blocks::ScoreBlock for #name #generics_ty #generics_where {
            fn column_schema(out: &mut crate::scoring::blocks::SchemaSink) {
                #(#schema_calls)*
            }

            fn columns(&self, out: &mut crate::scoring::blocks::ColSink) {
                #(#column_calls)*
            }

            fn linear_features(&self, out: &mut crate::scoring::blocks::FrameSink) {
                #(#linear_feature_calls)*
            }

            fn linear_feature_names(out: &mut crate::scoring::blocks::NameSink) {
                #(#linear_name_calls)*
            }

            fn nonlinear_features(&self, out: &mut crate::scoring::blocks::FrameSink) {
                #(#nonlinear_feature_calls)*
            }

            fn nonlinear_feature_names(out: &mut crate::scoring::blocks::NameSink) {
                #(#nonlinear_name_calls)*
            }
        }
    })
}

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
    fn emits_six_methods_and_routes_lane() {
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
        assert!(ts.contains("fn linear_features"));
        assert!(ts.contains("fn linear_feature_names"));
        assert!(ts.contains("fn nonlinear_features"));
        assert!(ts.contains("fn nonlinear_feature_names"));
        // a (log2) -> linear lane; b (raw) -> nonlinear lane; c -> column-only
        assert!(ts.contains("push_log2"));
        assert!(ts.contains("\"a\""));
        // all three fields appear in column_schema
        assert!(ts.contains("\"c\""));
    }

    #[test]
    fn rejects_unknown_generator() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(bogus)] pub a: f32 }
        })
        .unwrap();
        assert!(derive_score_block(di).is_err());
    }

    #[test]
    fn array_field_uses_slice_and_count() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(raw)] pub arr: [f32; NUM_MS2_IONS] }
        })
        .unwrap();
        let ts = derive_score_block(di).unwrap().to_string();
        assert!(ts.contains("push_slice"));
        assert!(ts.contains("f32_array"));
        assert!(ts.contains("NUM_MS2_IONS"));
    }

    #[test]
    fn rejects_non_raw_generator_on_array() {
        let di: syn::DeriveInput = syn::parse2(quote! {
            pub struct T { #[feat(log2)] pub arr: [f32; NUM_MS2_IONS] }
        })
        .unwrap();
        assert!(derive_score_block(di).is_err());
    }
}
