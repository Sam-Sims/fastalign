use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    author, version, about = "Quick multiple sequnce alignment using minimap2", long_about = None
)]
pub struct Cli {
    /// Input (unaligned) FASTA file.
    #[arg(
        short = 'i', long = "input", value_name = "Unaligned FASTA", required_unless_present = "help", value_parser(check_input_exists)
    )]
    pub input: String,

    /// Input reference FASTA file.
    #[arg(
        short = 'r', long = "reference", value_name = "Reference FASTA", required_unless_present = "help", value_parser(check_input_exists)
    )]
    pub reference: String,

    /// Output alignment file.
    #[arg(
        short = 'o', long = "output", value_name = "Output FASTA",
    )]
    pub output: String,

    /// Number of threads to use.
    /// Default: 1
    #[arg(short = 't', long = "threads", value_name = "Threads", default_value = "1")]
    pub threads: usize,

}

fn check_input_exists(s: &str) -> Result<String, String> {
    if s == "-" {
        return Ok(s.to_string());
    }
    if std::path::Path::new(s).exists() {
        Ok(s.to_string())
    } else {
        Err(format!("File does not exist: {}", s))
    }
}