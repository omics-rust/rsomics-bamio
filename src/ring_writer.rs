//! Parallel BGZF output with a fixed ring of reusable buffers.
//!
//! Compression workers reuse slots; a writer thread emits completed blocks in
//! input order.

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::num::NonZero;
use std::thread::JoinHandle;

use crossbeam_channel::{Receiver, Sender, bounded};
use libdeflater::{CompressionLvl, Compressor, Crc};
use rsomics_common::{Result, RsomicsError};

// The headroom keeps the complete BGZF member within the u16 BSIZE limit.
const MAX_IDATA: usize = 65_280 - 18 - 8 - 15; // = 65_239

const MAX_CDATA: usize = MAX_IDATA + 10;

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const HEADER_SIZE: usize = 18;
const TRAILER_SIZE: usize = 8;

struct Slot {
    idata: Vec<u8>,
    cdata: Vec<u8>,
    crc32: u32,
    clen: usize,
}

impl Slot {
    fn new() -> Self {
        Slot {
            idata: Vec::with_capacity(MAX_IDATA),
            cdata: vec![0u8; MAX_CDATA],
            crc32: 0,
            clen: 0,
        }
    }
}

type WorkItem = (usize, u64, usize);
type DoneItem = (usize, u64, usize, usize, u32);
type FreeItem = usize;

/// A deterministic parallel BGZF writer with reusable block buffers.
pub struct RingBgzfWriter<W: Write + Send + 'static> {
    slots: Vec<std::sync::Arc<std::sync::Mutex<Slot>>>,
    cur: usize,
    seq: u64,
    work_tx: Option<Sender<WorkItem>>,
    free_rx: Receiver<FreeItem>,
    done_tx: Option<Sender<DoneItem>>,
    deflate_handles: Vec<JoinHandle<()>>,
    writer_handle: Option<JoinHandle<io::Result<W>>>,
    finished: bool,
}

impl<W: Write + Send + 'static> RingBgzfWriter<W> {
    /// Create a writer with `workers` deflate threads over `sink`.
    pub fn new(sink: W, workers: NonZero<usize>) -> Self {
        let n_slots = workers.get() * 2 + 1;

        let slots: Vec<_> = (0..n_slots)
            .map(|_| std::sync::Arc::new(std::sync::Mutex::new(Slot::new())))
            .collect();

        let (work_tx, work_rx) = bounded::<WorkItem>(workers.get());
        let (free_tx, free_rx) = bounded::<FreeItem>(n_slots);
        for i in 1..n_slots {
            free_tx.send(i).unwrap();
        }

        let (done_tx, done_rx) = bounded::<DoneItem>(n_slots);

        let mut deflate_handles = Vec::with_capacity(workers.get());
        for _ in 0..workers.get() {
            let work_rx = work_rx.clone();
            let done_tx = done_tx.clone();
            let slots = slots.clone();
            deflate_handles.push(std::thread::spawn(move || {
                deflate_worker(slots, work_rx, done_tx);
            }));
        }

        let slots_w = slots.clone();
        let free_tx_w = free_tx;
        let writer_handle =
            std::thread::spawn(move || write_worker(sink, slots_w, done_rx, free_tx_w));

        RingBgzfWriter {
            slots,
            cur: 0,
            seq: 0,
            work_tx: Some(work_tx),
            free_rx,
            done_tx: Some(done_tx),
            deflate_handles,
            writer_handle: Some(writer_handle),
            finished: false,
        }
    }

    /// Finish the BGZF stream and return the underlying writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.flush_current()?;
        self.shutdown()
    }

    fn flush_current(&mut self) -> io::Result<()> {
        let idata_len = {
            let slot = self.slots[self.cur].lock().unwrap();
            slot.idata.len()
        };
        if idata_len == 0 {
            return Ok(());
        }
        let idx = self.cur;
        let seq = self.seq;
        self.seq += 1;
        self.work_tx
            .as_ref()
            .ok_or_else(|| io::Error::other("BGZF writer is finished"))?
            .send((idx, seq, idata_len))
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "deflate worker died"))?;

        let next = self
            .free_rx
            .recv()
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "writer thread died"))?;
        {
            let mut slot = self.slots[next].lock().unwrap();
            slot.idata.clear();
        }
        self.cur = next;
        Ok(())
    }

    fn shutdown(&mut self) -> io::Result<W> {
        if self.finished {
            return Err(io::Error::other("BGZF writer is already finished"));
        }
        self.finished = true;

        drop(self.work_tx.take());
        let mut worker_panicked = false;
        for handle in self.deflate_handles.drain(..) {
            worker_panicked |= handle.join().is_err();
        }
        drop(self.done_tx.take());

        let writer = self
            .writer_handle
            .take()
            .ok_or_else(|| io::Error::other("BGZF writer thread is missing"))?
            .join()
            .map_err(|_| io::Error::other("BGZF writer thread panicked"))??;

        if worker_panicked {
            Err(io::Error::other("BGZF compression worker panicked"))
        } else {
            Ok(writer)
        }
    }
}

impl<W: Write + Send + 'static> Write for RingBgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let remaining = {
            let slot = self.slots[self.cur].lock().unwrap();
            MAX_IDATA - slot.idata.len()
        };
        let amt = remaining.min(buf.len());
        {
            let mut slot = self.slots[self.cur].lock().unwrap();
            slot.idata.extend_from_slice(&buf[..amt]);
        }
        let is_full = {
            let slot = self.slots[self.cur].lock().unwrap();
            slot.idata.len() >= MAX_IDATA
        };
        if is_full {
            self.flush_current()?;
        }
        Ok(amt)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_current()
    }
}

impl<W: Write + Send + 'static> Drop for RingBgzfWriter<W> {
    fn drop(&mut self) {
        if !self.finished {
            let _ = self.flush_current();
            let _ = self.shutdown();
        }
    }
}

fn deflate_worker(
    slots: Vec<std::sync::Arc<std::sync::Mutex<Slot>>>,
    work_rx: Receiver<WorkItem>,
    done_tx: Sender<DoneItem>,
) {
    let mut compressor = Compressor::new(CompressionLvl::new(6).expect("level 6 is valid"));
    let mut idata_scratch: Vec<u8> = Vec::with_capacity(MAX_IDATA);

    while let Ok((idx, seq, idata_len)) = work_rx.recv() {
        let (clen, crc32) = {
            let mut slot = slots[idx].lock().unwrap();

            idata_scratch.resize(idata_len, 0);
            idata_scratch.copy_from_slice(&slot.idata[..idata_len]);

            let mut crc_engine = Crc::new();
            crc_engine.update(&idata_scratch);
            let crc32 = crc_engine.sum();

            let clen = compressor
                .deflate_compress(&idata_scratch, &mut slot.cdata)
                .expect("deflate compress cannot fail on valid input");
            slot.clen = clen;
            slot.crc32 = crc32;
            (clen, crc32)
        };
        if done_tx.send((idx, seq, idata_len, clen, crc32)).is_err() {
            break;
        }
    }
}

fn write_worker<W: Write + Send + 'static>(
    mut sink: W,
    slots: Vec<std::sync::Arc<std::sync::Mutex<Slot>>>,
    done_rx: Receiver<DoneItem>,
    free_tx: Sender<FreeItem>,
) -> io::Result<W> {
    let mut pending: BTreeMap<u64, (usize, usize, usize, u32)> = BTreeMap::new();
    let mut next_seq: u64 = 0;
    let mut frame_buf: Vec<u8> = Vec::with_capacity(HEADER_SIZE + MAX_CDATA + TRAILER_SIZE);

    while let Ok((idx, seq, idata_len, clen, crc32)) = done_rx.recv() {
        pending.insert(seq, (idx, idata_len, clen, crc32));

        while let Some(&(pidx, pidata_len, pclen, pcrc32)) = pending.get(&next_seq) {
            pending.remove(&next_seq);
            next_seq += 1;

            {
                let slot = slots[pidx].lock().unwrap();
                write_bgzf_block(
                    &mut sink,
                    &mut frame_buf,
                    &slot.cdata[..pclen],
                    pcrc32,
                    pidata_len,
                )?;
            }
            if free_tx.send(pidx).is_err() {
                break;
            }
        }
    }

    while let Some((&seq, _)) = pending.first_key_value() {
        if seq != next_seq {
            break;
        }
        let (pidx, pidata_len, pclen, pcrc32) = pending.remove(&next_seq).unwrap();
        next_seq += 1;
        let slot = slots[pidx].lock().unwrap();
        write_bgzf_block(
            &mut sink,
            &mut frame_buf,
            &slot.cdata[..pclen],
            pcrc32,
            pidata_len,
        )?;
    }

    sink.write_all(&BGZF_EOF)?;
    sink.flush()?;
    Ok(sink)
}

fn write_bgzf_block<W: Write>(
    sink: &mut W,
    frame_buf: &mut Vec<u8>,
    cdata: &[u8],
    crc32: u32,
    idata_len: usize,
) -> io::Result<()> {
    let block_size = HEADER_SIZE + cdata.len() + TRAILER_SIZE;

    frame_buf.clear();

    frame_buf.extend_from_slice(&[
        0x1f, 0x8b, // ID1, ID2
        0x08, // CM = DEFLATE
        0x04, // FLG = FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME = 0
        0x00, // XFL = 0
        0xff, // OS = 255 (unknown)
        0x06, 0x00, // XLEN = 6
        b'B', b'C', // SI1, SI2
        0x02, 0x00, // SLEN = 2
    ]);
    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
    frame_buf.extend_from_slice(&bsize.to_le_bytes());

    frame_buf.extend_from_slice(cdata);
    frame_buf.extend_from_slice(&crc32.to_le_bytes());
    let isize_u32 = u32::try_from(idata_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
    frame_buf.extend_from_slice(&isize_u32.to_le_bytes());

    sink.write_all(frame_buf)
}

/// Create a [`RingBgzfWriter`] over a new file.
pub fn create_ring_bgzf<P: AsRef<std::path::Path>>(
    path: P,
    workers: NonZero<usize>,
) -> Result<RingBgzfWriter<std::fs::File>> {
    let file = std::fs::File::create(path.as_ref()).map_err(|e| {
        RsomicsError::InvalidInput(format!("creating {}: {e}", path.as_ref().display()))
    })?;
    Ok(RingBgzfWriter::new(file, workers))
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};
    use std::sync::{Arc, Mutex};

    use super::*;

    #[derive(Clone)]
    struct SharedBuf(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedBuf {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn new_buf() -> (SharedBuf, Arc<Mutex<Vec<u8>>>) {
        let inner = Arc::new(Mutex::new(Vec::new()));
        (SharedBuf(Arc::clone(&inner)), inner)
    }

    #[test]
    fn round_trip_single_thread() {
        let workers = NonZero::new(1).unwrap();
        let (sink, data) = new_buf();
        let mut w = RingBgzfWriter::new(sink, workers);
        w.write_all(b"hello world").unwrap();
        w.finish().unwrap();

        let bytes = data.lock().unwrap().clone();
        let mut reader = noodles::bgzf::io::Reader::new(std::io::Cursor::new(bytes));
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).unwrap();
        assert_eq!(decompressed, b"hello world");
    }

    #[test]
    fn round_trip_multi_thread() {
        for nw in [2usize, 4, 8] {
            let workers = NonZero::new(nw).unwrap();
            let payload: Vec<u8> = (0u8..=255).cycle().take(MAX_IDATA * 4 + 1000).collect();
            let (sink, data) = new_buf();
            let mut w = RingBgzfWriter::new(sink, workers);
            w.write_all(&payload).unwrap();
            w.finish().unwrap();

            let bytes = data.lock().unwrap().clone();
            let mut reader = noodles::bgzf::io::Reader::new(std::io::Cursor::new(bytes));
            let mut decompressed = Vec::new();
            reader.read_to_end(&mut decompressed).unwrap();
            assert_eq!(decompressed, payload, "failed at workers={nw}");
        }
    }

    #[test]
    fn eof_block_present() {
        let workers = NonZero::new(2).unwrap();
        let (sink, data) = new_buf();
        let mut w = RingBgzfWriter::new(sink, workers);
        w.write_all(b"test").unwrap();
        w.finish().unwrap();

        let bytes = data.lock().unwrap().clone();
        let end = bytes.len();
        assert_eq!(&bytes[end - BGZF_EOF.len()..], BGZF_EOF);
    }

    #[test]
    fn output_identical_across_thread_counts() {
        let payload: Vec<u8> = (0u8..=255).cycle().take(MAX_IDATA * 3 + 500).collect();

        let outputs: Vec<Vec<u8>> = [1usize, 2, 4]
            .iter()
            .map(|&nw| {
                let workers = NonZero::new(nw).unwrap();
                let (sink, data) = new_buf();
                let mut w = RingBgzfWriter::new(sink, workers);
                w.write_all(&payload).unwrap();
                w.finish().unwrap();
                data.lock().unwrap().clone()
            })
            .collect();

        assert_eq!(outputs[0], outputs[1], "t1 vs t2 differ");
        assert_eq!(outputs[0], outputs[2], "t1 vs t4 differ");
    }

    #[test]
    fn finish_surfaces_sink_flush_error() {
        #[derive(Debug)]
        struct FlushError;

        impl Write for FlushError {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                Ok(buf.len())
            }

            fn flush(&mut self) -> io::Result<()> {
                Err(io::Error::other("flush failed"))
            }
        }

        let mut writer = RingBgzfWriter::new(FlushError, NonZero::<usize>::MIN);
        writer.write_all(b"test").unwrap();
        let error = writer.finish().unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::Other);
        assert_eq!(error.to_string(), "flush failed");
    }
}
