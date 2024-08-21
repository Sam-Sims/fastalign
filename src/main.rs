use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::str::FromStr;
use clap::Parser;
use minimap2::*;
use noodles::fasta;
use anyhow::{Result, Context, anyhow};
use noodles::fasta::Record;
use noodles::fasta::record::{Definition, Sequence};
use std::thread;
use crossbeam_channel::{unbounded};

mod cli;

enum CigarOperation {
    Match(usize),
    Insertion(usize),
    Deletion(usize),
    Skipped(usize),
    SoftClip(usize),
    HardClip(usize),
    Pad(usize),
    Equal(usize),
    Diff(usize),
}

/// FromStr implementation for CigarOperation parsing
impl FromStr for CigarOperation {
    type Err = anyhow::Error;

    fn from_str(cigar_operation: &str) -> Result<Self, Self::Err> {
        // split the string into the count and operation
        // count = number of times to apply the operation
        let count = cigar_operation[..cigar_operation.len()-1].parse::<usize>()
            .context(format!("Failed to parse CIGAR operation count: {}", cigar_operation))?;
        // operation = the type of operation to apply
        match cigar_operation.chars().last().context(format!("Failed to parse CIGAR operation character: {}", cigar_operation))? {
            'M' => Ok(CigarOperation::Match(count)),
            'I' => Ok(CigarOperation::Insertion(count)),
            'D' => Ok(CigarOperation::Deletion(count)),
            'N' => Ok(CigarOperation::Skipped(count)),
            'S' => Ok(CigarOperation::SoftClip(count)),
            'H' => Ok(CigarOperation::HardClip(count)),
            'P' => Ok(CigarOperation::Pad(count)),
            '=' => Ok(CigarOperation::Equal(count)),
            'X' => Ok(CigarOperation::Diff(count)),
            c => Err(anyhow!("Unknown CIGAR operation: {} in CIGAR string {}", c, cigar_operation)),
        }
        // notes from SAM spec:
        // H can only be present as the first and/or last operation. TODO: Sanity check this?
        // For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments, the interpretation of N is not defined. TODO: Error on N?
        // Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ. TODO: Sanity check this?
    }
}

/// Split the cigar into individual operations and parse
fn parse_cigar(cigar_string: &str) -> Result<Vec<CigarOperation>> {
    cigar_string.split_inclusive(char::is_alphabetic)
        .map(|cigar_operation| cigar_operation.parse().with_context(|| format!("Failed to parse CIGAR operation: {}", cigar_operation)))
        .collect()
}

/// Build an aligned sequence from the CIGAR string
fn align_sequence(sequence: &[u8], reference_len: usize, cigar: &str, aln_start: i32) -> Result<Vec<u8>> {
    let mut aligned_seq = Vec::with_capacity(reference_len);
    // Add gaps for any reference bases before the start of the alignment
    aligned_seq.extend("-".repeat(aln_start as usize).bytes());

    let mut seq_pos = 0;
    let mut ref_pos = aln_start as usize;

    for op in parse_cigar(cigar)? {
        // Process the CIAGAR operations
        // Currently only handles M, I, D, and N operations
        // Insertions are ignored to ensure each sequence matches the ref length
        // TODO: Handle insertions and output modified reference
        match op {
            CigarOperation::Match(count) | CigarOperation::Equal(count) | CigarOperation::Diff(count) => {
                let end_pos = seq_pos + count;
                if end_pos > sequence.len() {
                    return Err(anyhow!(
                        "CIGAR operation out-of-bounds sequence: seq_pos={}, count={}, sequence length={}",
                        seq_pos, count, sequence.len()
                    ));
                }
                //aligned_seq.push_str(&sequence[seq_pos..end_pos]);
                aligned_seq.extend_from_slice(&sequence[seq_pos..end_pos]);
                seq_pos += count;
                ref_pos += count;
            },
            // Insertions are ignored, so just increment the sequence position
            CigarOperation::Insertion(count) => seq_pos += count,
            // Might get Ns in CIGAR?
            CigarOperation::Deletion(count) | CigarOperation::Skipped(count) => {
                aligned_seq.extend("-".repeat(count).bytes());
                ref_pos += count;
            },
            // TODO: Soft and hard clips are ignored for now, decide how to handle them
            CigarOperation::SoftClip(count) | CigarOperation::HardClip(count) => seq_pos += count,
            _ => {}
        }
    }

    // Add gaps for any reference bases after the end of the alignment
    if ref_pos < reference_len {
        aligned_seq.extend("-".repeat(reference_len - ref_pos).bytes());
    }

    Ok(aligned_seq)
}


fn align_record(record: &fasta::Record, reference: &str, aligner: &Aligner) -> Result<fasta::Record> {
    let seq = record.sequence();
    let name = record.name();

    let alignment = aligner.map(seq.as_ref(), false, false, None, None)
        .map_err(|e| anyhow!(e))
        .context("Failed to align sequence")?;
    
    // 

    // alignment should always contain only 1 alignment, but if mapping fails, it might be empty
    if let Some(aln) = alignment.first() {
        if !aln.is_primary {
            println!("Not a primary alignment: {}", std::str::from_utf8(name)?);
        }
        if let Some(cigar_string) = aln.alignment.as_ref().and_then(|a| a.cigar_str.as_ref()) {
            let aligned_seq = align_sequence(
                seq.as_ref(),
                reference.len(),
                cigar_string,
                aln.target_start,
            ).context("Failed to align sequence")?;

            let definition = Definition::new(name.to_owned(), None);
            let sequence = Sequence::from(aligned_seq);
            Ok(Record::new(definition, sequence))
        } else {
            Err(anyhow!("No CIGAR string found for alignment {}", std::str::from_utf8(name)?))
        }
    } else {
        Err(anyhow!("No alignment found for sequence {}", std::str::from_utf8(name)?))
    }
}


fn process_fasta(input_path: &str, output_path: &str, reference: &str, aligner: &Aligner, num_threads: usize) -> Result<()> {
    let input_file = File::open(input_path).context("Failed to open input file")?;
    let mut input_reader = fasta::Reader::new(BufReader::new(input_file));

    let output_file = File::create(output_path).context("Failed to create output file")?;
    let mut output_writer = fasta::Writer::new(BufWriter::new(output_file));

    let (record_snd, record_recv) = unbounded();
    let (aligned_snd, aligned_recv) = unbounded();

    thread::scope(|s| -> Result<()> {
        // Spawn a thread to read the input FASTA file and send records to record_snd
        s.spawn(|| -> Result<()> {
            for record in input_reader.records() {
                let record = record.context("Failed to read FASTA record")?;
                record_snd.send(record).context("Failed to send record")?;
            }
            drop(record_snd);
            Ok(())
        });

        // Create threads to receive records from record_recv, align them, and send them to aligned_snd
        for _ in 0..num_threads {
            let record_receiver = record_recv.clone();
            let result_sender = aligned_snd.clone();
            let mut aligner = aligner.clone();

            s.spawn(move || -> Result<()> {
                while let Ok(record) = record_receiver.recv() {
                    let aligned_record = align_record(&record, &reference, &aligner)
                        .context("Failed to align record")?;
                    result_sender.send(aligned_record).context("Failed to send aligned record")?;
                }
                // stops a double free seg fault, see https://github.com/jguhlin/minimap2-rs/issues/71
                aligner.idx = None;
                Ok(())
            });
        }
        drop(aligned_snd);

        // Final thread to receive aligned records from aligned_recv and write them
        s.spawn(move || -> Result<()>{
            while let Ok(aligned_record) = aligned_recv.recv() {
                output_writer.write_record(&aligned_record)
                    .context("Failed to write aligned record")?;
            }
            Ok(())
        });

        Ok(())
    }).context("Thread error")?;

    Ok(())
}

fn fastalign() -> Result<()> {
    let args = cli::Cli::parse();

    let aligner = Aligner {
        mapopt: MapOpt {
            sc_ambi: 0,
            ..Aligner::builder().asm20().mapopt
        },
        ..Aligner::builder().asm20()
    }
        .with_cigar()
        .with_sam_hit_only()
        .with_index(&args.reference, None)
        .map_err(|e| anyhow!(e))
        .context("Failed to build aligner")?;

    let mut reference = String::new();
    let ref_file = File::open(&args.reference).context("Failed to open reference file")?;
    let mut ref_reader = fasta::Reader::new(BufReader::new(ref_file));
    for record in ref_reader.records() {
        let record = record.context("Failed to read reference FASTA record")?;
        reference.push_str(std::str::from_utf8(record.sequence().as_ref()).context("Invalid UTF-8 reference sequence")?);
    }

    let num_threads = args.threads;
    process_fasta(&args.input, &args.output, &reference, &aligner, num_threads)?;

    Ok(())
}

fn main() {
    if let Err(e) = fastalign() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}