use super::poly::Polynomial;


pub struct PolyMatrix<const K: usize, const L: usize> {
    data: [[Polynomial; L]; K]
}

impl <const K: usize, const L: usize>PolyMatrix<K, L> {
    pub fn zeroes() -> Self {
        let data = [[Polynomial::default(); L]; K];
        Self {data}
    }


    pub fn set(&mut self, i: usize, j: usize, poly: Polynomial) {
        let stored_poly = self.data.get_mut(i).expect("").get_mut(j).expect("");
        *stored_poly = poly;
    }
}

impl <const K: usize, const L: usize>Default for PolyMatrix<K, L> {
    fn default() -> Self {
        Self::zeroes()
    }
}