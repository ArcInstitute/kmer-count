use anyhow::{bail, Result};
use clap::Parser;

#[derive(Parser)]
pub struct Cli {
    /// The input file to process (*.{fq.gz, fq, bq, vbq})
    #[clap(required = true)]
    input: String,

    /// The k-mer size to count
    #[clap(short = 'k', long, default_value = "4")]
    kmer: usize,

    /// The number of threads to use
    ///
    /// 0 for maximum available
    #[clap(short = 'T', long, default_value = "1")]
    threads: usize,
}
impl Cli {
    pub fn format(&self) -> Result<FileFormat> {
        let ext = if self.input.contains(".gz") {
            let num_splits = self.input.split(".").count();
            self.input.split(".").nth(num_splits - 2)
        } else {
            self.input.split(".").last()
        };
        if let Some(ext) = ext {
            match ext {
                "bq" => Ok(FileFormat::Bq),
                "vbq" => Ok(FileFormat::Vbq),
                "fq" | "fastq" => Ok(FileFormat::Fq),
                _ => bail!(
                    "Could not determine file format from {} with ext {}",
                    &self.input,
                    ext
                ),
            }
        } else {
            bail!("Could not determine file format from {}", &self.input);
        }
    }
    pub fn path(&self) -> &str {
        &self.input
    }
    pub fn ksize(&self) -> usize {
        self.kmer
    }
    pub fn threads(&self) -> usize {
        if self.threads == 0 {
            num_cpus::get()
        } else {
            self.threads.min(num_cpus::get())
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
pub enum FileFormat {
    Bq,
    Vbq,
    Fq,
}
