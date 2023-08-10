use super::poly::Polynomial;

#[derive(Debug)]
pub struct PolyMatrix<const K: usize, const L: usize> {
    data: [[Polynomial; L]; K]
}

impl <const K: usize, const L: usize>PolyMatrix<K, L> {
    pub fn zeroes() -> Self {
        let data = [[Polynomial::default(); L]; K];
        Self {data}
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut Polynomial {
        let poly = self.data.get_mut(i).unwrap().get_mut(j).unwrap();
        poly
    }
}

impl <const K: usize, const L: usize>Default for PolyMatrix<K, L> {
    fn default() -> Self {
        Self::zeroes()
    }
}