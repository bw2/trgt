use std::cmp::Ordering;

#[inline]
fn median_of_three_index(data: &[i32], low: usize, mid: usize, high: usize) -> usize {
    let a = data[low];
    let b = data[mid];
    let c = data[high];
    if (a <= b && b <= c) || (c <= b && b <= a) {
        mid // b is the median
    } else if (b <= a && a <= c) || (c <= a && a <= b) {
        low // a is the median
    } else {
        high // c is the median
    }
}

/// Partition data in-place using Lomuto scheme
fn partition_inplace(data: &mut [i32], low: usize, high: usize) -> usize {
    if low >= high {
        return low;
    }
    let mid = low + (high - low) / 2;
    let pivot_index = median_of_three_index(data, low, mid, high);

    // Swap the chosen pivot element to the end of the range [low..=high]
    data.swap(pivot_index, high);
    let pivot_value = data[high];

    // `i` tracks the boundary between elements <= pivot and elements > pivot
    let mut i = low;
    for j in low..high {
        if data[j] <= pivot_value {
            data.swap(i, j);
            i += 1;
        }
    }
    data.swap(i, high);
    i
}

// Iterative Quickselect algorithm using median-of-three to avoid worst case quadratic runtime
fn select_inplace(data: &mut [i32], k: usize) -> Option<i32> {
    if data.is_empty() || k >= data.len() {
        return None;
    }

    let mut low = 0;
    let mut high = data.len() - 1;

    loop {
        if low == high {
            return if low == k { Some(data[low]) } else { None };
        }
        let pivot_index = partition_inplace(data, low, high);
        match pivot_index.cmp(&k) {
            Ordering::Equal => return Some(data[k]),
            Ordering::Greater => {
                if pivot_index == 0 {
                    return None;
                }
                high = pivot_index - 1;
            }
            Ordering::Less => {
                low = pivot_index + 1;
            }
        }
        if low > high {
            return None;
        }
    }
}

pub fn median(data: &[i32]) -> Option<f32> {
    let size = data.len();
    if size == 0 {
        return None;
    }
    let mut data_copy = data.to_vec();
    match size {
        even if even % 2 == 0 => {
            let k1 = (even / 2) - 1;
            let k2 = even / 2;
            let fst_med_opt = select_inplace(&mut data_copy, k1);
            match fst_med_opt {
                Some(fst) => {
                    let min_in_right_partition = data_copy[k2..].iter().min();
                    min_in_right_partition.map(|&snd| (fst + snd) as f32 / 2.0)
                }
                None => None,
            }
        }
        odd => {
            let k = odd / 2;
            select_inplace(&mut data_copy, k).map(|x| x as f32)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rng, seq::SliceRandom};

    fn naive_median(data: &[i32]) -> Option<f32> {
        let size = data.len();
        if size == 0 {
            return None;
        }
        let mut sorted_data = data.to_vec();
        sorted_data.sort();
        if size % 2 == 0 {
            let mid1 = sorted_data[(size / 2) - 1];
            let mid2 = sorted_data[size / 2];
            Some((mid1 + mid2) as f32 / 2.0)
        } else {
            Some(sorted_data[size / 2] as f32)
        }
    }

    #[test]
    fn test_median_empty() {
        let data: [i32; 0] = [];
        assert_eq!(median(&data), None);
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_single_element() {
        let data = [5];
        assert_eq!(median(&data), Some(5.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_odd_count_unsorted() {
        let data = [3, 1, 4, 1, 5];
        assert_eq!(median(&data), Some(3.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_even_count_unsorted() {
        let data = [3, 1, 4, 2];
        assert_eq!(median(&data), Some(2.5));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_odd_count_sorted() {
        let data = [1, 2, 3, 4, 5];
        assert_eq!(median(&data), Some(3.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_even_count_sorted() {
        let data = [1, 2, 3, 4];
        assert_eq!(median(&data), Some(2.5));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_odd_count_reverse_sorted() {
        let data = [5, 4, 3, 2, 1];
        assert_eq!(median(&data), Some(3.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_even_count_reverse_sorted() {
        let data = [4, 3, 2, 1];
        assert_eq!(median(&data), Some(2.5));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_with_duplicates_even() {
        let data = [1, 2, 2, 3];
        assert_eq!(median(&data), Some(2.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_with_duplicates_odd() {
        let data = [1, 2, 2, 2, 3];
        assert_eq!(median(&data), Some(2.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_with_negatives_even() {
        let data = [-1, -5, 0, 2];
        assert_eq!(median(&data), Some(-0.5));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_with_negatives_odd() {
        let data = [-1, -5, 0, 2, -3];
        assert_eq!(median(&data), Some(-1.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_all_same() {
        let data = [7, 7, 7, 7, 7];
        assert_eq!(median(&data), Some(7.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_large_even() {
        let data: Vec<i32> = (1..=1000).collect();
        assert_eq!(median(&data), Some(500.5));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_large_odd() {
        let data: Vec<i32> = (1..=999).collect();
        assert_eq!(median(&data), Some(500.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_large_reverse_sorted_even() {
        let data: Vec<i32> = (1..=1000).rev().collect();
        assert_eq!(median(&data), Some(500.5));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_large_reverse_sorted_odd() {
        let data: Vec<i32> = (1..=999).rev().collect();
        assert_eq!(median(&data), Some(500.0));
        assert_eq!(median(&data), naive_median(&data));
    }

    #[test]
    fn test_median_large_random() {
        let mut data: Vec<i32> = (1..=1001).collect();
        data.shuffle(&mut rng());
        assert_eq!(median(&data), Some(501.0));
        assert_eq!(median(&data), naive_median(&data));

        let mut data_even: Vec<i32> = (1..=1000).collect();
        data_even.shuffle(&mut rng());
        assert_eq!(median(&data_even), Some(500.5));
        assert_eq!(median(&data_even), naive_median(&data_even));
    }

    #[test]
    fn test_select_inplace_basic() {
        let data = vec![3, 1, 4, 1, 5, 9, 2, 6];
        assert_eq!(select_inplace(&mut data.clone(), 0), Some(1));
        assert_eq!(select_inplace(&mut data.clone(), 7), Some(9));
        assert_eq!(select_inplace(&mut data.clone(), 3), Some(3));
        assert_eq!(select_inplace(&mut data.clone(), 4), Some(4));
    }

    #[test]
    fn test_select_inplace_duplicates() {
        let data = vec![3, 1, 4, 1, 5, 1, 9];
        assert_eq!(select_inplace(&mut data.clone(), 0), Some(1));
        assert_eq!(select_inplace(&mut data.clone(), 1), Some(1));
        assert_eq!(select_inplace(&mut data.clone(), 2), Some(1));
        assert_eq!(select_inplace(&mut data.clone(), 3), Some(3));
    }

    #[test]
    fn test_select_inplace_sorted() {
        let data = vec![1, 2, 3, 4, 5, 6];
        assert_eq!(select_inplace(&mut data.clone(), 2), Some(3));
        assert_eq!(select_inplace(&mut data.clone(), 5), Some(6));
    }

    #[test]
    fn test_select_inplace_reverse_sorted() {
        let data = vec![6, 5, 4, 3, 2, 1];
        assert_eq!(select_inplace(&mut data.clone(), 0), Some(1));
        assert_eq!(select_inplace(&mut data.clone(), 3), Some(4));
    }

    #[test]
    fn test_select_inplace_out_of_bounds() {
        let data = vec![3, 1, 4];
        assert_eq!(select_inplace(&mut data.clone(), 3), None);
        assert_eq!(select_inplace(&mut data.clone(), 10), None);
        let mut empty_data: Vec<i32> = vec![];
        assert_eq!(select_inplace(&mut empty_data, 0), None);
    }
}
