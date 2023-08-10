use crate::algebra::poly::Polynomial;

pub struct PolyVec<const N: usize> {
    data: [Polynomial; N]
}

impl <const N: usize>Default for PolyVec<N> {
    fn default() -> Self {
        let data = [Polynomial::default(); N];
        Self {
            data
        }
    }
}

impl <const N: usize>IntoIterator for PolyVec<N> {
    type Item = Polynomial;
    type IntoIter = std::array::IntoIter<Self::Item, N>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}