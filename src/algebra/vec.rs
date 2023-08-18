use crate::algebra::poly::Polynomial;
use std::ops::{AddAssign, Index, IndexMut, Mul, Sub, SubAssign};
use std::slice::IterMut;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PolyVec<const N: usize> {
    pub(crate) data: [Polynomial; N],
}

impl<const N: usize> PolyVec<N> {
    pub fn iter_mut(&mut self) -> IterMut<Polynomial> {
        self.data.iter_mut()
    }

    pub fn to_ntt(&self) -> Self {
        let mut clone = self.clone();
        for i in 0..N {
            clone[i].ntt();
        }
        clone
    }

    pub fn caddq(&mut self) {
        for i in 0..N {
            self[i].caddq();
        }
    }

    pub fn reduce(&mut self) {
        for i in 0..N {
            self[i].reduce()
        }
    }

    /// This function returns a Boolean indicating whether the infinite norm of the vector of
    /// polynomials is greater than or equal to a specific bound.
    ///
    /// # Arguments:
    /// * `bound` - The upper bound
    pub fn check_infinite_norm(&self, bound: usize) -> bool {
        for poly in &self.data {
            if poly.check_infinite_norm(bound) {
                return true;
            }
        }

        false
    }

    pub fn power2round(&self) -> (PolyVec<N>, PolyVec<N>) {
        let mut poly_vec_zero = PolyVec::<N>::default();
        let mut poly_vec_one = PolyVec::<N>::default();

        for i in 0..N {
            let (poly_zero, poly_one) = self[i].power2round();
            poly_vec_zero[i] = poly_zero;
            poly_vec_one[i] = poly_one;
        }

        (poly_vec_zero, poly_vec_one)
    }

    pub fn expand_mask<const ARR_SIZE: usize, const GAMMA_ONE: usize>(
        rho_prime: &[u8],
        kappa: u16,
    ) -> PolyVec<N> {
        let mut vec = PolyVec::default();
        let mut nonce = kappa;

        for i in 0..N {
            Polynomial::fill_uniform_gamma::<ARR_SIZE, GAMMA_ONE>(rho_prime, nonce, &mut vec[i]);
            nonce += 1;
        }

        vec
    }

    pub fn decompose<const GAMMA_TWO: usize>(vec: &PolyVec<N>) -> (PolyVec<N>, PolyVec<N>) {
        let mut r_zero = PolyVec::<N>::default();
        let mut r_one = PolyVec::<N>::default();

        for i in 0..N {
            Polynomial::decompose::<GAMMA_TWO>(&mut r_zero[i], &mut r_one[i], &vec[i]);
        }

        return (r_one, r_zero);
    }

    /// This function fills a hint vector based on low part of an input vector.
    ///
    /// # Arguments
    /// * `output_vec` - The hint vector to be filled
    /// * `v_zero` - The low part of an input vector
    /// * `v_one` - The high part of an input vector
    pub fn make_hint<const GAMMA_TWO: usize>(
        output_vec: &mut PolyVec<N>,
        v_zero: &PolyVec<N>,
        v_one: &PolyVec<N>,
    ) -> usize {
        let mut number_ones = 0;

        for i in 0..N {
            number_ones +=
                Polynomial::make_hint::<GAMMA_TWO>(&mut output_vec[i], &v_zero[i], &v_one[i]);
        }

        number_ones
    }

    pub fn expand_s<const ETA: usize>(rho_prime: &[u8], nonce_init_value: u16) -> PolyVec<N> {
        let mut vec = PolyVec::<N>::default();
        let mut nonce = nonce_init_value;
        vec.iter_mut().for_each(|poly| {
            Polynomial::fill_uniform_eta::<ETA>(rho_prime, nonce, poly);
            nonce += 1
        });
        vec
    }

    /// This function multiplies the vector of polynomials of length K by 2^D without modular
    /// reduction.
    ///
    /// The function assumes input coefficients to be less that 2^{31 - D}.
    pub fn shiftl(&mut self) {
        for i in 0..N {
            self[i].shiftl();
        }
    }

    /// This function uses the hint vector to correct the high bits of input vector and outputs
    /// a vector filled with the correction.
    ///
    /// # Arguments
    ///
    /// * `h` - A hint vector
    /// * `vec` - A vector of polynomials for which we should correct the high bits.
    pub fn use_hint<const K: usize, const GAMMA_TWO: usize>(
        h: &PolyVec<K>,
        vec: &PolyVec<K>,
    ) -> PolyVec<K> {
        let mut w_one_prime = PolyVec::<K>::default();

        for i in 0..K {
            Polynomial::use_hint::<GAMMA_TWO>(&h[i], &vec[i], &mut w_one_prime[i])
        }

        w_one_prime
    }
}

impl<const N: usize> Default for PolyVec<N> {
    fn default() -> Self {
        let data = [Polynomial::default(); N];
        Self { data }
    }
}

impl<const N: usize> Index<usize> for PolyVec<N> {
    type Output = Polynomial;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<const N: usize> IndexMut<usize> for PolyVec<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<const N: usize> AddAssign<&PolyVec<N>> for PolyVec<N> {
    fn add_assign(&mut self, rhs: &PolyVec<N>) {
        for i in 0..N {
            self[i] += &rhs[i]
        }
    }
}

impl<const N: usize> Mul<&PolyVec<N>> for PolyVec<N> {
    type Output = Polynomial;
    fn mul(self, rhs: &PolyVec<N>) -> Self::Output {
        let mut poly = Polynomial::default();

        for i in 0..N {
            poly += &(self[i] * &rhs[i]);
        }

        poly.reduce();
        poly.ntt_to_montgomery();
        poly
    }
}

impl<const N: usize> Mul<&Polynomial> for PolyVec<N> {
    type Output = PolyVec<N>;

    fn mul(self, rhs: &Polynomial) -> Self::Output {
        let mut vec = PolyVec::<N>::default();

        for i in 0..N {
            vec[i] = self[i] * rhs;
            vec[i].ntt_to_montgomery();
            vec[i].reduce();
        }

        vec
    }
}

impl<const N: usize> Sub<&PolyVec<N>> for PolyVec<N> {
    type Output = PolyVec<N>;
    fn sub(self, rhs: &PolyVec<N>) -> Self::Output {
        let mut vec = PolyVec::default();

        for i in 0..N {
            vec[i] = self[i] - &rhs[i]
        }

        vec
    }
}

impl<const N: usize> SubAssign<&PolyVec<N>> for PolyVec<N> {
    fn sub_assign(&mut self, rhs: &PolyVec<N>) {
        for i in 0..N {
            self[i] = self[i] - &rhs[i]
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::algebra::poly::Polynomial;
    use crate::algebra::vec::PolyVec;
    use crate::constants::{H_TEST, Q, W_ONE_TEST, W_TEST, W_ZERO_TEST};

    #[test]
    fn test_should_produce_the_correct_h_vec() {
        let mut h = PolyVec::<4>::default();
        const LOW_PART_VEC: PolyVec<4> = PolyVec {
            data: [
                Polynomial {
                    coefficients: [
                        89961, -67247, -72071, -20390, 13181, 33288, 45560, -95780, -29401, 82826,
                        -12376, 45452, 38779, -80617, -9829, 27912, -107600, -50088, -98914, 45568,
                        -39289, 26306, 13373, 44384, 23112, 80731, -40655, 25809, -2928, -9888,
                        38102, -57881, -4726, -91950, 55795, -124385, 43373, -70733, -42841, 46946,
                        -34519, -4023, 35974, -89931, -81864, 104098, 15084, -68123, 60726, -22959,
                        37156, -27470, 18821, 111439, -50187, -24378, 28953, 20333, -15257, -80768,
                        18308, 1740, -44628, -91626, -94703, 31987, 83184, 62844, 71889, -48340,
                        43658, 4470, 57459, -38768, 19377, -65666, 35134, -69613, -35283, 38055,
                        -35880, 45787, 18101, -67891, -92165, 84346, 15201, -59982, -104917, 41490,
                        -17077, -112615, -9738, -97666, 44587, 69967, -35809, -76140, 24420,
                        -36581, 8795, -76887, -79841, 77629, -5418, -91486, 53909, 22878, -67299,
                        -13878, 97868, 51400, 38403, -41979, -35804, -75473, -29945, 47354, 5837,
                        8469, -27172, -100116, -19213, 29748, 16764, 48293, -19506, -22085, 79512,
                        -88388, 42997, -7472, -80470, 104666, 85313, 3662, -38290, 75970, 82419,
                        6997, 36062, 31583, -28848, 46064, 49954, -57, 100696, -46722, 1575,
                        -43176, 16871, 28084, 60246, -41358, -21691, -81263, -15989, -92696,
                        -35366, 18910, 33045, 19555, 65489, 29023, -44552, 48339, 4093, 55042,
                        -87561, 89171, 87363, 25288, 38104, -119824, -61387, -34567, -48313, 88241,
                        15468, 14254, 26259, 43593, -36664, -101112, 95660, 24578, 56496, 12531,
                        55834, 104450, -8105, -69022, -59961, 2057, -57015, 87267, -7739, 84711,
                        10987, -26934, -57535, 50852, 12728, 12420, 15429, -78383, 88585, 104509,
                        40304, 76301, 48445, -74399, -31954, 51402, 13589, -11316, -31026, 13293,
                        -56819, 63607, 109322, 24282, 43052, 73273, 64363, -46336, 3044, 116883,
                        -45288, 70726, -33638, -28508, -25213, 39053, -62819, 71100, 61529, 7617,
                        -47351, -75259, 48159, 71960, 25290, -34869, 56344, -4238, 76071, 14881,
                        -48102, 106932, 70126, 44835, 49886, 14162, 43572, -70633,
                    ],
                },
                Polynomial {
                    coefficients: [
                        -8097, -29640, -29680, -33890, -15518, 4693, -7860, -66577, -18819, 20679,
                        -98313, 17526, -25780, -35326, 52529, -46309, 87481, -13657, 68596, 22959,
                        -30477, 773, -42711, -95738, -3125, 22087, -37808, -58657, -21981, 7545,
                        -35869, 102677, -26083, -69801, -62109, 48911, 74083, -14301, 98402, 46443,
                        -23239, 58834, -103041, -69352, -91949, 16323, 90747, -55016, 58083,
                        -14534, -42385, -69283, -29519, 61411, -65093, -57062, -95321, -35681,
                        -12624, -65530, -69451, -64674, 87073, -94429, 27666, 13238, 77783, 82044,
                        -100309, -69707, 27050, -50631, 71668, 24259, 25389, 28822, 790, -92773,
                        -49196, 48565, -74865, 47503, -22736, 904, 72710, -73909, -26745, 73275,
                        -76285, -80642, -31490, -26293, -25846, 47034, -87270, 85170, 91393, 15974,
                        -99502, 73232, 67287, -70698, -56599, 49991, -10539, 69838, 5787, -17257,
                        54359, 75291, -64378, -88783, 70444, -97343, 114038, 21448, 9483, -21118,
                        -46188, 67886, -69138, 90902, -27247, 64926, -92542, 45486, 31657, 2881,
                        25955, -103473, 76982, -25701, -23224, 77797, -13840, -70509, 98407,
                        -55156, 56925, -7584, -51202, 60652, 33471, 13827, -75985, 72955, 14831,
                        -9879, -88201, 56423, 27649, -46233, -52761, -7103, 24588, -89049, -47294,
                        -30260, 69258, 85202, 11760, 48098, -69568, -4645, -75774, -8524, 70479,
                        49287, -67384, -24231, -75576, -41778, -50595, 22433, 43932, 25458, -18342,
                        8392, -125400, -48734, -14895, 63506, 9011, -32953, -37706, -20169, 2950,
                        -54343, -15305, 24298, 81296, 20374, 83472, 18708, 50931, -50702, -53534,
                        -21187, 20762, -52398, -68244, 49473, -76071, -38126, -80385, 13797,
                        -60673, -27195, 54337, -69589, 45012, -71585, 53912, 55700, -9676, 21316,
                        -102318, 2872, -64948, 68773, -20003, -30451, -88318, -74136, 50539,
                        -11406, 14165, 46063, 4990, 41982, -30706, 18414, 100380, -62623, 38978,
                        -77861, -64009, 80513, -98884, -54617, 51361, -75938, -44567, 83407,
                        -53813, 82192, -33500, 35698, 32055, -50042, 24966, 8349, -31102, 68077,
                        -78846, -20160,
                    ],
                },
                Polynomial {
                    coefficients: [
                        24601, -47298, -95785, 59656, -55750, 51777, 9529, -21998, -8914, -20445,
                        26573, 87179, -58954, -64443, 15925, 40577, 44201, 252, 64905, -20023,
                        -5352, 20823, 73873, -19201, 48504, 25413, 9104, 66932, 16599, 62045,
                        75934, -10049, -5291, 66524, -19354, 75425, -4236, -43668, -69726, -112849,
                        -28926, -25135, -86104, 26928, -137985, -9379, 30662, -104851, 67409, 9514,
                        -155, 64637, -34472, -14657, 29370, 59471, 7377, -54537, 46583, 81705,
                        -56599, 42202, 56920, 18973, 110270, 16020, -5448, -68009, -30282, -45755,
                        32823, 479, -41238, 14695, -5059, 20533, -82399, -111212, -33548, 1220,
                        -5375, -84791, -81557, 33606, -45214, -68960, -5923, 3757, 23842, -28343,
                        42840, -19593, -56207, -38165, 58472, -47264, -38327, 60554, -71299,
                        -56107, -13907, -2312, -19498, 112669, -51115, 47545, -28840, -25374,
                        -94718, -95206, -46747, -90882, -38076, -2753, -21093, -3591, 8954, 36990,
                        82703, -48745, -3532, -83139, 93119, 81223, 84006, -22290, -79097, 37775,
                        -97942, 60557, 91745, 86903, 45137, 96413, 69179, -90773, -69109, 8137,
                        88878, -8202, -77106, -60238, -83264, -72403, -52360, 5357, -64026, 66168,
                        -20861, -11894, -89949, -4542, 55675, 72174, 50056, 46333, 50629, -85987,
                        36955, 50389, 23192, -94551, -57933, 40869, -48866, 48342, -60176, -5937,
                        -42196, 71328, -77626, 70692, -37103, -52282, -91398, 79788, 33927, 34506,
                        104988, 101080, -28872, -86060, 52754, 85951, 55201, -58672, -56853, 51551,
                        -107525, -68262, -117608, 8288, 49078, 83007, -16759, -106176, 64379,
                        20742, 7993, -43377, -4074, -69031, 4285, -15307, -56630, 57958, -17367,
                        11832, -13430, 23260, 37093, -95029, -92736, -10029, -92213, -37480,
                        -43777, 12847, 4856, 17153, 109186, 45503, 23273, -2960, 46829, 84964,
                        -7974, -15012, 67727, 77730, 97923, 16222, -3176, 12954, -93474, -66539,
                        -81241, 11620, 33150, -81266, 79205, -107592, -51936, 130013, 32112, 12229,
                        -79179, 36257, 65744, -25700, -90241, 82002, -16148, -42441, 81025, 53709,
                    ],
                },
                Polynomial {
                    coefficients: [
                        9706, -55515, 86676, -14556, -77355, 92334, 69282, 69546, -4752, 54512,
                        84299, -30325, 57393, 58277, -87948, -96666, -11098, 12074, 33388, -10897,
                        4632, -64626, 88563, -7732, -48897, -27449, 35096, 51675, 16627, 61709,
                        55599, -24291, 37840, 95034, -43214, 10758, 22768, -5774, 27358, -9444,
                        15184, -54787, -18121, 79991, -45663, 33597, 2540, 12455, -1606, 80480,
                        62312, 40987, -2924, 17460, -22656, 91284, 1446, 57780, 32869, 28501,
                        75253, 1110, 66988, -61848, 70920, -61337, -43814, -32882, -74188, 5153,
                        -20446, 56505, -88954, -87115, 9729, 6000, 19021, 60743, 94854, 7395,
                        49085, -28842, 45346, 90471, 26185, -84444, -10990, -106665, 46826, 66827,
                        -81746, 79243, -22984, 28341, -79814, 3821, -81964, -81728, -100529,
                        -91296, -68962, 73367, -82758, 68387, -27143, 82730, 84871, 56550, 28923,
                        3664, 38851, 78845, 84807, 76768, -103846, 64844, 81366, 36934, -42369,
                        -44509, 13984, 1119, 50821, -10644, -35212, 15249, -57114, -6979, 90729,
                        79034, 101485, 92642, -69166, -94532, 73599, 75320, -96381, -28551, 54347,
                        -59800, 93815, 62766, -28100, 46475, -48952, -40121, 615, 55294, -9338,
                        66483, 108, -43195, 40755, 25935, -91366, 30590, -41141, 70184, -89995,
                        84623, -9405, 1252, -65667, 47772, -10076, 23525, -4999, 50286, 5136,
                        -64262, -97197, 56430, 2669, -3213, 32815, -32861, -18815, 80512, 64939,
                        -38766, 65282, -15960, -73548, -88701, -98523, -21835, 64018, 83724, 36207,
                        -88616, 74732, -64203, 28368, -80644, 12840, 79932, 44753, -46839, -89212,
                        -15772, -54368, -21148, -57999, -52038, -61884, 15079, 19374, -29111,
                        74231, 12608, 84075, 29937, 18423, -22994, -67708, 100127, -45177, 63764,
                        -76722, -59376, -18659, -4462, -39890, 65375, -88532, -45665, 45446,
                        -72112, -91432, -2912, -74496, 47634, -22188, 37071, -9242, 93795, -612,
                        5360, 63159, -89291, 89337, 76599, 2232, -57554, 74022, 33109, -50136,
                        -100615, 28642, -15053, 19883, 34708, 15209, 94840, -107182, -34613,
                    ],
                },
            ],
        };

        PolyVec::<4>::make_hint::<{ (Q as usize - 1) / 88 }>(&mut h, &LOW_PART_VEC, &W_ONE_TEST);

        println!("H: {:?}", h);
        assert_eq!(h[0].coefficients, H_TEST[0].coefficients);
        assert_eq!(h[1].coefficients, H_TEST[1].coefficients);
        assert_eq!(h[2].coefficients, H_TEST[2].coefficients);
        assert_eq!(h[3].coefficients, H_TEST[3].coefficients);
    }

    #[test]
    fn test_should_decompose_correctly_the_w_vector() {
        let (w_one, w_zero) = PolyVec::<4>::decompose::<{ (Q as usize - 1) / 88 }>(&W_TEST);

        assert_eq!(w_one, W_ONE_TEST);
        assert_eq!(w_zero, W_ZERO_TEST);
    }
}
