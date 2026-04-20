use crate::utils::math::lnfact_f32;

#[inline(always)]
pub fn single_lazyscore(slc: impl IntoIterator<Item = f32>) -> f32 {
    lnfact_f32(slc.into_iter().map(|x| x.max(1.0).ln()).sum())
}
