use std::{
    fs::File,
    io::{stdout, BufWriter, Write},
};

use anyhow::Result;
use clap::Parser;

mod cli;
mod count;
use cli::{Cli, FileFormat};
use count::Counter;
use paraseq::{fastq, parallel::ParallelReader};
use serde::Serialize;

#[derive(Serialize)]
pub struct Statistics {
    /// Number of elapsed seconds
    elapsed: f64,
    /// Number of records processed
    n_records: usize,
    /// File format processed
    format: FileFormat,
    /// Number of threads used
    threads: usize,
    /// Throughput
    throughput: f64,
    /// Input path
    path: String,
}
impl Statistics {
    pub fn new(
        elapsed: f64,
        n_records: usize,
        format: FileFormat,
        threads: usize,
        path: String,
    ) -> Self {
        Self {
            elapsed,
            n_records,
            format,
            threads,
            path,
            throughput: n_records as f64 / elapsed,
        }
    }
}

fn main() -> Result<()> {
    let args = Cli::parse();
    let mut out_handle = BufWriter::new(stdout());
    let start = std::time::Instant::now();
    let counter = match args.format()? {
        FileFormat::Bq => {
            let reader = binseq::MmapReader::new(args.path())?;
            let counter = Counter::new(args.ksize())?;
            reader.process_parallel(counter.clone(), args.threads())?;
            counter
        }
        FileFormat::Vbq => {
            let reader = vbinseq::MmapReader::new(args.path())?;
            let counter = Counter::new(args.ksize())?;
            reader.process_parallel(counter.clone(), args.threads())?;
            counter
        }
        FileFormat::Fq => {
            let (in_handle, _comp) = niffler::send::from_path(args.path())?;
            let reader = fastq::Reader::new(in_handle);
            let counter = Counter::new(args.ksize())?;
            reader.process_parallel(counter.clone(), args.threads())?;
            counter
        }
    };
    if let Some(log_path) = args.log_path() {
        let elapsed = start.elapsed().as_secs_f64();
        let statistics = Statistics::new(
            elapsed,
            counter.n_records(),
            args.format()?,
            args.threads(),
            args.path().to_string(),
        );
        let mut log_handle = File::create(log_path)?;
        serde_json::to_writer_pretty(&mut log_handle, &statistics)?;
        log_handle.flush()?;
    }

    counter.pprint(&mut out_handle)?;

    Ok(())
}
