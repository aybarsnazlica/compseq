use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation::*;
use bio::alignment::pairwise::*;
use bio::scores::blosum62;
use bio::utils::TextSlice;
use std::cmp::max;

const GAP_OPEN_PENALTY: i32 = -5;
const GAP_EXTEND_PENALTY: i32 = -1;

pub fn align(x: TextSlice, y: TextSlice, mode: &str) -> Alignment {
    let mut aligner = Aligner::with_capacity(
        x.len(),
        y.len(),
        GAP_OPEN_PENALTY,
        GAP_EXTEND_PENALTY,
        &blosum62,
    );

    let alignment = match mode {
        "global" => aligner.global(x, y),
        "local" => aligner.local(x, y),
        _ => panic!("Select either global or local alignment"),
    };
    alignment
}

pub fn identity(alignment: &Alignment) -> u8 {
    let num_matches = alignment
        .operations
        .iter()
        .filter(|&&op| op == Match)
        .count();

    let length = alignment.xlen.max(alignment.ylen);

    (num_matches as f64 / length as f64 * 100.0).round() as u8
}

pub fn similarity(alignment: &Alignment, x: TextSlice, y: TextSlice) -> u8 {
    let (query, reference) = if x.len() > y.len() { (y, x) } else { (x, y) };
    let mut score = 0;
    let mut max_score = 0;

    // Left unaligned region
    max_score += reference[..alignment.ystart]
        .iter()
        .map(|&c| blosum62(c, c))
        .sum::<i32>();

    let mut xi = alignment.xstart;
    let mut yi = alignment.ystart;

    for op in &alignment.operations {
        match op {
            Match => {
                if xi < query.len() && yi < y.len() {
                    let s = blosum62(query[xi], reference[yi]);
                    score += s;
                    max_score += s;
                }

                xi += 1;
                yi += 1;
            }
            Subst => {
                if xi < query.len() && yi < y.len() {
                    score += blosum62(query[xi], reference[yi]);
                    max_score += max(
                        blosum62(query[xi], query[xi]),
                        blosum62(reference[yi], reference[yi]),
                    );
                }

                xi += 1;
                yi += 1;
            }
            Del => {
                score += GAP_OPEN_PENALTY;
                max_score += blosum62(query[xi], query[xi]);
                xi += 1;
            }
            Ins => {
                score += GAP_OPEN_PENALTY;
                max_score += blosum62(reference[yi], reference[yi]);
                yi += 1;
            }
            Xclip(n) => xi += n,
            Yclip(n) => yi += n,
        }
    }

    // Right unaligned region
    max_score += reference[alignment.yend..]
        .iter()
        .map(|&c| blosum62(c, c))
        .sum::<i32>();

    if max_score == 0 {
        return 0;
    }

    ((score.max(0) as f64 / max_score as f64) * 100.0).round() as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_operations() {
        let x = b"LSPADKTNVKAA";
        let y = b"PEEKSAV";
        let alignment = align(x, y, "global");
        assert_eq!(
            alignment.operations,
            [
                Ins, Ins, Match, Subst, Subst, Match, Subst, Ins, Ins, Ins, Match, Subst
            ]
        );
    }

    #[test]
    fn identity_exact_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"LSPADKTNVKAA";
        let alignment = align(x, y, "global");
        assert_eq!(identity(&alignment), 100);
    }

    #[test]
    fn identity_no_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"QQQQQQQQQQQQ";
        let alignment = align(x, y, "global");
        assert_eq!(identity(&alignment), 0);
    }

    #[test]
    fn identity_ninety() {
        let x = b"LSPADKTNVK";
        let y = b"LSPADQTNVK";
        let alignment = align(x, y, "global");
        assert_eq!(identity(&alignment), 90);
    }

    #[test]
    fn identity_zero() {
        let x = b"LSPADKTNVK";
        let y = b"";
        let alignment = align(x, y, "global");
        assert_eq!(identity(&alignment), 0);
    }

    #[test]
    fn similarity_exact_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"LSPADKTNVKAA";
        let alignment = align(x, y, "global");
        assert_eq!(similarity(&alignment, x, y), 100);
    }

    #[test]
    fn similarity_no_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"QQQQQQQQQQ";
        let alignment = align(x, y, "global");
        assert_eq!(similarity(&alignment, x, y), 0);
    }

    #[test]
    fn similarity_single_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"LLLLL";
        let alignment = align(x, y, "global");
        assert_eq!(similarity(&alignment, x, y), 0);
    }

    #[test]
    fn similarity_zero() {
        let x = b"LSPADKTNVKAA";
        let y = b"";
        let alignment = align(x, y, "global");
        assert_eq!(similarity(&alignment, x, y), 0);
    }

    #[test]
    fn similarity_ninety_two() {
        let x = b"LSPADKTNVK";
        let y = b"LSPADQTNVK";
        let alignment = align(x, y, "global");
        assert_eq!(similarity(&alignment, x, y), 92);
    }

    #[test]
    fn similarity_seventy_four() {
        let x = b"ALSPADQTNVK";
        let y = b"LSPADQTNVK";
        let alignment = align(x, y, "global");
        assert_eq!(similarity(&alignment, x, y), 74);
    }
}
