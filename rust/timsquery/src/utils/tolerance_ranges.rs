use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IncludedRange<T: Copy + PartialOrd>(T, T);

impl<T> IncludedRange<T>
where
    T: Copy + PartialOrd,
{
    pub fn new(left: T, right: T) -> Self {
        // Swap if left > right
        if left > right {
            Self(right, left)
        } else {
            Self(left, right)
        }
    }

    // TODO: Implement other overlaps ... (union ...)
    pub fn intersection(&self, other: Self) -> Option<Self> {
        let left = match self.start().partial_cmp(&other.start()) {
            Some(std::cmp::Ordering::Equal) => self.start(),
            Some(std::cmp::Ordering::Less) => other.start(),
            Some(std::cmp::Ordering::Greater) => self.start(),
            None => self.start(),
        };
        let right = match self.end().partial_cmp(&other.end()) {
            Some(std::cmp::Ordering::Equal) => self.start(),
            Some(std::cmp::Ordering::Less) => self.start(),
            Some(std::cmp::Ordering::Greater) => other.start(),
            None => self.start(),
        };

        if left > right {
            None
        } else {
            Some(Self(left, right))
        }
    }

    pub fn contains(&self, x: T) -> bool {
        self.0 <= x && x <= self.1
    }

    pub fn intersects(&self, x: Self) -> bool {
        self.contains(x.start()) || self.contains(x.end())
    }

    pub fn start(&self) -> T {
        self.0
    }

    pub fn end(&self) -> T {
        self.1
    }
}

impl<T: Copy + PartialOrd> From<(T, T)> for IncludedRange<T> {
    fn from(x: (T, T)) -> Self {
        Self::new(x.0, x.1)
    }
}

impl<T: Copy + PartialOrd> From<IncludedRange<T>> for (T, T) {
    fn from(val: IncludedRange<T>) -> Self {
        (val.0, val.1)
    }
}

pub fn ppm_tol_range(elem: f64, tol_ppm: f64) -> IncludedRange<f64> {
    let utol = elem * (tol_ppm / 1e6);
    let left_e = elem - utol;
    let right_e = elem + utol;
    (left_e, right_e).into()
}

pub fn pct_tol_range(elem: f64, tol_pct: f64) -> IncludedRange<f64> {
    let utol = elem * (tol_pct / 100.0);
    let left_e = elem - utol;
    let right_e = elem + utol;
    (left_e, right_e).into()
}

pub fn tof_tol_range(tof: u32, tol_ppm: f64, converter: &Tof2MzConverter) -> IncludedRange<u32> {
    let mz = converter.convert(tof);
    let mz_range = ppm_tol_range(mz, tol_ppm);
    IncludedRange::new(
        converter.invert(mz_range.start()).round() as u32,
        converter.invert(mz_range.end()).round() as u32,
    )
}

pub fn scan_tol_range(scan: u16, tol_pct: f64, converter: &Scan2ImConverter) -> IncludedRange<u16> {
    let im = converter.convert(scan as f64);
    let im_range = pct_tol_range(im, tol_pct);
    let scan_min = converter.invert(im_range.start()).round() as u16;
    let scan_max = converter.invert(im_range.end()).round() as u16;
    // Note I need to do this here bc the conversion between scan numbers and ion
    // mobilities is not monotonically increasing. IN OTHER WORDS, lower scan numbers
    // are higher 1/k0.... But im not sure if they are ALWAYS inversely proportional.
    if scan_min > scan_max {
        (scan_max, scan_min).into()
    } else {
        (scan_min, scan_max).into()
    }
}
