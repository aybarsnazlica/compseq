use compseq;

use clap::{Parser, Subcommand};

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

            let sequences = compseq::read_fasta(input);

            for (i, rec1) in sequences.iter().enumerate() {
                for rec2 in sequences.iter().skip(i + 1) {
                    let x = rec1.seq();
                    let y = rec2.seq();
                    let alignment = compseq::align(x, y);

                    println!(
                        "Alignment between {} and {}: score {}, identity {}",
                        rec1.id(),
                        rec2.id(),
                        alignment.score,
                        compseq::identity(&alignment)
                    );
                }
            }
        }
    }
}
