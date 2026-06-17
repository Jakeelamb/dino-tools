use std::io;

use dino_seq::{
    FastqChunkConfig, FastqChunkSinkExt, FastqReader, FastqRecordSink, FastqVisitRecord, Result,
};

#[derive(Debug)]
struct OwnedRead {
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

#[derive(Debug, Default)]
struct ReadBuffer {
    reads: Vec<OwnedRead>,
    bases: u64,
}

impl FastqRecordSink for ReadBuffer {
    fn record(&mut self, record: FastqVisitRecord<'_>) -> Result<()> {
        self.bases += record.seq().len() as u64;
        self.reads.push(OwnedRead {
            name: record.name().to_vec(),
            seq: record.seq().to_vec(),
            qual: record.qual().to_vec(),
        });
        Ok(())
    }
}

impl FastqChunkSinkExt for ReadBuffer {
    fn reserve_records(&mut self, records: usize) {
        self.reads.reserve(records);
    }
}

impl ReadBuffer {
    fn clear(&mut self) {
        self.reads.clear();
        self.bases = 0;
    }
}

fn main() -> Result<()> {
    let stdin = io::stdin();
    let mut reader = FastqReader::new(stdin.lock());
    let config = FastqChunkConfig::new(10_000_000).min_records(40_000);
    let mut sink = ReadBuffer::default();

    sink.reserve_records(config.estimated_records(150));
    while let Some(stats) = reader.next_chunk_with_sink(config, &mut sink)? {
        let checksum: usize = sink
            .reads
            .iter()
            .map(|read| read.name.len() ^ read.seq.len() ^ read.qual.len())
            .sum();
        eprintln!(
            "chunk first_record={} records={} bases={} checksum={}",
            stats.first_record_index(),
            stats.records(),
            stats.bases(),
            checksum
        );

        // Downstream tools typically hand `sink.reads` to alignment or another
        // processing stage here, then reuse the allocation for the next chunk.
        sink.clear();
        sink.reserve_records(config.estimated_records(150));
    }

    Ok(())
}
