use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation::*;
use bio::alignment::pairwise::*;
use bio::io::fasta;
use bio::scores::blosum62;
use bio::utils::TextSlice;
use std::cmp::max;
use std::fs;

const GAP_OPEN_PENALTY: i32 = -5;
const GAP_EXTEND_PENALTY: i32 = -1;

pub fn align(x: TextSlice, y: TextSlice) -> Alignment {
    let mut aligner = Aligner::with_capacity(
        x.len(),
        y.len(),
        GAP_OPEN_PENALTY,
        GAP_EXTEND_PENALTY,
        &blosum62,
    );
    let alignment = aligner.global(x, y);
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
    let mut score = 0;
    let mut max_score = 0;
    let mut xi = alignment.xstart;
    let mut yi = alignment.ystart;

    for op in &alignment.operations {
        match op {
            Match | Subst => {
                if xi < alignment.xend && yi < alignment.yend {
                    score += blosum62(x[xi], y[yi]);
                    let score_aa = blosum62(x[xi], x[xi]);
                    let score_bb = blosum62(y[yi], y[yi]);
                    max_score += max(score_aa, score_bb);
                }
                xi += 1;
                yi += 1;
            }
            Del => {
                if xi < alignment.xend {
                    score += GAP_OPEN_PENALTY;
                    max_score += blosum62(x[xi], x[xi]);
                }
                xi += 1;
            }
            Ins => {
                if yi < alignment.yend {
                    score += GAP_OPEN_PENALTY;
                    max_score += blosum62(y[yi], y[yi]);
                }
                yi += 1;
            }
            Xclip(n) => xi += n,
            Yclip(n) => yi += n,
        }
    }

    if max_score == 0 {
        return 0;
    }

    ((score.max(0) as f64 / max_score as f64) * 100.0).round() as u8
}

pub fn read_fasta(fasta_file: String) -> Vec<fasta::Record> {
    let file = fs::File::open(fasta_file).expect("File opening failed");
    let records = fasta::Reader::new(file).records();

    let sequences: Vec<fasta::Record> = records.map(|r| r.expect("Invalid FASTA record")).collect();
    sequences
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
            [
                Ins, Ins, Match, Subst, Subst, Match, Subst, Ins, Ins, Ins, Match, Subst
            ]
        );
    }

    #[test]
    fn identity_exact_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"LSPADKTNVKAA";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 100);
    }

    #[test]
    fn identity_no_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"QQQQQQQQQQQQ";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 0);
    }

    #[test]
    fn identity_ninety() {
        let x = b"LSPADKTNVK";
        let y = b"LSPADQTNVK";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 90);
    }

    #[test]
    fn identity_zero() {
        let x = b"LSPADKTNVK";
        let y = b"";
        let alignment = align(x, y);
        assert_eq!(identity(&alignment), 0);
    }

    #[test]
    fn similarity_exact_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"LSPADKTNVKAA";
        let alignment = align(x, y);
        assert_eq!(similarity(&alignment, x, y), 100);
    }

    #[test]
    fn similarity_no_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"QQQQQQQQQQ";
        let alignment = align(x, y);
        assert_eq!(similarity(&alignment, x, y), 0);
    }

    #[test]
    fn similarity_single_match() {
        let x = b"LSPADKTNVKAA";
        let y = b"LLLLL";
        let alignment = align(x, y);
        assert_eq!(similarity(&alignment, x, y), 0);
    }

    #[test]
    fn similarity_zero() {
        let x = b"LSPADKTNVKAA";
        let y = b"";
        let alignment = align(x, y);
        assert_eq!(similarity(&alignment, x, y), 0);
    }

    #[test]
    fn similarity_ninety_two() {
        let x = b"LSPADKTNVK";
        let y = b"LSPADQTNVK";
        let alignment = align(x, y);
        assert_eq!(similarity(&alignment, x, y), 92);
    }

    // TODO: Handle this case
    // #[test]
    // fn similarity_ninety_two() {
    //     let x = b"ALSPADKTNVK";
    //     let y = b"LSPADQTNVK";
    //     let alignment = align(x, y);
    //     for ali in &alignment.operations {
    //         println!("{:#?}", ali);
    //     }
    //     println!("{}", alignment.pretty(x, y, 80));
    //     assert_eq!(similarity(&alignment, x, y), 92);
    // }
}
