use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation::*;
use bio::alignment::pairwise::*;
use bio::scores::blosum62;
use bio::utils::TextSlice;

pub fn align(x: TextSlice, y: TextSlice) -> Alignment {
    let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &blosum62);
    let alignment = aligner.local(x, y);
    alignment
}

pub fn identity(alignment: &Alignment) -> f64 {
    let num_matches = alignment
        .operations
        .iter()
        .filter(|&&op| op == Match)
        .count() as f64;

    let length = alignment.xlen.max(alignment.ylen) as f64;

    num_matches / length
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_operations() {
        let x = b"LSPADKTNVKAA";
        let y = b"PEEKSAV";
        let alignment = align(x, y);
        assert_eq!(
            alignment.operations,
            [Match, Subst, Subst, Match, Subst, Subst, Match]
        );
    }

    #[test]
    fn identity_exact_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"LSPADKTNVKAA";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 1.0);
    }

    #[test]
    fn identity_no_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"QQQQQQQQQQQQ";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 0.0);
    }

    #[test]
    fn identity_ninety() {
        let x = b"LSPADKTNVK";
        let y = b"LSPADQTNVK";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 0.9);
    }

    #[test]
    fn identity_zero() {
        let x = b"LSPADKTNVK";
        let y = b"";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 0.0);
    }
}
