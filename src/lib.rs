use bio::utils::TextSlice;
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::scores::blosum62;


pub fn align(x: TextSlice, y: TextSlice) -> Alignment {
    let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &blosum62);
    let alignment = aligner.local(x, y);
    alignment
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::alignment::AlignmentOperation::*;

    #[test]
    fn alignment() {
        let x = b"LSPADKTNVKAA";
        let y = b"PEEKSAV";
        let alignment = align(x, y);
        assert_eq!(
    alignment.operations,
    [Match, Subst, Subst, Match, Subst, Subst, Match]
);
    }
}