use clap::{Parser, Subcommand};

use compseq::alignment::{align, identity, similarity};
use compseq::io::read_fasta;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Align {
        #[arg(short, long)]
        mode: String,

        #[arg(short, long)]
        input: String,
    },
}

fn main() {
    let args = Args::parse();

    match args.command {
        Command::Align { mode, input } => {
            println!("Running alignment in mode: {}", mode);

            let sequences = read_fasta(input);

            for (i, rec1) in sequences.iter().enumerate() {
                for rec2 in sequences.iter().skip(i + 1) {
                    let x = rec1.seq();
                    let y = rec2.seq();
                    let alignment = align(x, y, &mode);

                    println!(
                        "Alignment between {} and {}: identity {}%, similarity {}%.",
                        rec1.id(),
                        rec2.id(),
                        identity(&alignment),
                        similarity(&alignment, x, y)
                    );
                }
            }
        }
    }
}
