use std::{io::Write, sync::Arc};

use anyhow::{bail, Result};
use hashbrown::HashMap;
use parking_lot::Mutex;

pub type KmerMap = HashMap<u64, usize>;

#[derive(Clone)]
pub struct Counter {
    /// Maps
    local_map: KmerMap,
    global_map: Arc<Mutex<KmerMap>>,

    /// Decoding buffer (binseq)
    dbuf: Vec<u8>,

    /// Kmer-size
    ksize: usize,
}
impl Counter {
    pub fn new(ksize: usize) -> Result<Self> {
        if ksize > 32 {
            bail!("Currently only supporting kmer sizes <32bp.")
        }
        Ok(Self {
            local_map: KmerMap::default(),
            global_map: Arc::new(Mutex::new(KmerMap::default())),
            dbuf: Vec::new(),
            ksize,
        })
    }

    pub fn pprint<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut dbuf = Vec::with_capacity(self.ksize);
        for (k, v) in self.global_map.lock().iter() {
            // Decode the sequence
            {
                dbuf.clear();
                bitnuc::from_2bit(*k, self.ksize, &mut dbuf)?;
            }
            let seq_str = std::str::from_utf8(&dbuf)?;
            writeln!(writer, "{}\t{}", v, seq_str)?;
        }
        writer.flush()?;
        Ok(())
    }
}

fn kmer_count(seq: &[u8], map: &mut KmerMap, ksize: usize) -> Result<()> {
    if seq.len() < ksize {
        bail!(
            "Could not process record of size {} with ksize {}",
            seq.len(),
            ksize
        )
    }
    for idx in 0..seq.len() - ksize {
        let sub = &seq[idx..idx + ksize];
        if let Ok(encoding) = bitnuc::as_2bit(sub) {
            *map.entry(encoding).or_insert(0) += 1;
        } else {
            // Skip invalid nucleotides
        }
    }
    Ok(())
}

/// Takes all keys and entries from one kmer map and updates the other.
fn consolidate_maps(from: &mut KmerMap, to: &mut KmerMap) {
    from.drain()
        .for_each(|(k, v)| *to.entry(k).or_default() += v);
}

impl binseq::ParallelProcessor for Counter {
    fn process_record(&mut self, record: binseq::RefRecord) -> binseq::Result<()> {
        // Decode record
        {
            self.dbuf.clear();
            record.decode_s(&mut self.dbuf)?;
        }
        kmer_count(&self.dbuf, &mut self.local_map, self.ksize)?;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        // Lock the global map
        {
            let mut global = self.global_map.lock();
            consolidate_maps(&mut self.local_map, &mut global);
        }
        self.local_map.clear(); // should be drained but just in case
        Ok(())
    }
}

impl vbinseq::ParallelProcessor for Counter {
    fn process_record(&mut self, record: vbinseq::RefRecord) -> vbinseq::Result<()> {
        // Decode record
        {
            self.dbuf.clear();
            record.decode_s(&mut self.dbuf)?;
        }
        kmer_count(&self.dbuf, &mut self.local_map, self.ksize)?;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> vbinseq::Result<()> {
        // Lock the global map
        {
            let mut global = self.global_map.lock();
            consolidate_maps(&mut self.local_map, &mut global);
        }
        self.local_map.clear(); // should be drained but just in case
        Ok(())
    }
}

impl paraseq::parallel::ParallelProcessor for Counter {
    fn process_record<Rf: paraseq::fastx::Record>(
        &mut self,
        record: Rf,
    ) -> paraseq::parallel::Result<()> {
        kmer_count(record.seq(), &mut self.local_map, self.ksize)?;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Lock the global map
        {
            let mut global = self.global_map.lock();
            consolidate_maps(&mut self.local_map, &mut global);
        }
        self.local_map.clear(); // should be drained but just in case
        Ok(())
    }
}
