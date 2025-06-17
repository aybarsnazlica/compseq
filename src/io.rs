use bio::io::fasta;
use std::fs;

pub fn read_fasta(fasta_file: String) -> Vec<fasta::Record> {
    let file = fs::File::open(fasta_file).expect("File opening failed");
    let records = fasta::Reader::new(file).records();

    let sequences: Vec<fasta::Record> = records.map(|r| r.expect("Invalid FASTA record")).collect();
    sequences
}
