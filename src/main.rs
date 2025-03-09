use std::io::{stdout, BufWriter};

use anyhow::Result;
use clap::Parser;

mod cli;
mod count;
use cli::{Cli, FileFormat};
use count::Counter;
use paraseq::{fastq, parallel::ParallelReader};

fn main() -> Result<()> {
    let args = Cli::parse();
    let mut out_handle = BufWriter::new(stdout());
    match args.format()? {
        FileFormat::Bq => {
            let reader = binseq::MmapReader::new(args.path())?;
            let counter = Counter::new(args.ksize())?;
            reader.process_parallel(counter.clone(), args.threads())?;
            counter.pprint(&mut out_handle)
        }
        FileFormat::Vbq => {
            let reader = vbinseq::MmapReader::new(args.path())?;
            let counter = Counter::new(args.ksize())?;
            reader.process_parallel(counter.clone(), args.threads())?;
            counter.pprint(&mut out_handle)
        }
        FileFormat::Fq => {
            let (in_handle, _comp) = niffler::send::from_path(args.path())?;
            let reader = fastq::Reader::new(in_handle);
            let counter = Counter::new(args.ksize())?;
            reader.process_parallel(counter.clone(), args.threads())?;
            counter.pprint(&mut out_handle)
        }
    }?;
    Ok(())
}
