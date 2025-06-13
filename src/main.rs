use clap::{Parser, Subcommand};
use compseq::align;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Align two sequences
    Align {
        #[arg(short, long)]
        mode: String,
    },
}

fn main() {
    let args = Args::parse();

    match args.command {
        Command::Align { mode } => {
            println!("Running alignment in mode: {}", mode);
            
            let x = b"LSPADKTNVKAA";
            let y = b"PEEKSAV";
            let alignment = align(x, y);
            println!("Alignment score {}", alignment.score)
        }
    }
}